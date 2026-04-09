from __future__ import annotations

import argparse
import json
import math
import re
from pathlib import Path
from typing import Iterable

import pandas as pd


ROOT = Path(__file__).resolve().parent
LIMMA_DIR = ROOT / "limma_confirmation_output"
OUTPUT_DIR = ROOT / "ortholog_enrichment_output"

INTERACTION_INPUT = LIMMA_DIR / "interaction_effects_limma.csv"
FLOWERING_INPUT = LIMMA_DIR / "flowering_candidates_limma.csv"

INTERACTION_MAP = OUTPUT_DIR / "interaction_gene_ortholog_map.csv"
FLOWERING_MAP = OUTPUT_DIR / "flowering_interaction_ortholog_map.csv"
ALL_GPROFILER = OUTPUT_DIR / "all_interactions_gprofiler.csv"
FLOWERING_GPROFILER = OUTPUT_DIR / "flowering_interactions_gprofiler.csv"
ALL_PLOT = OUTPUT_DIR / "all_interactions_enrichment.png"
FLOWERING_PLOT = OUTPUT_DIR / "flowering_interactions_enrichment.png"
SUMMARY_JSON = OUTPUT_DIR / "ortholog_enrichment_summary.json"
REPORT_MD = OUTPUT_DIR / "ortholog_enrichment_report.md"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Recover and rerun ortholog-based enrichment analysis for limma "
            "interaction genes."
        )
    )
    parser.add_argument(
        "--reuse-existing-maps",
        action="store_true",
        help="Reuse saved ortholog map CSVs if they already exist.",
    )
    parser.add_argument(
        "--skip-remap",
        action="store_true",
        help="Do not call MyGene. Requires existing ortholog map CSVs.",
    )
    parser.add_argument(
        "--top-n",
        type=int,
        default=12,
        help="Number of enriched terms to show in plots and markdown tables.",
    )
    return parser.parse_args()


def clean_text(value: object) -> str:
    if pd.isna(value):
        return ""
    return str(value).strip()


def looks_like_locus(value: str) -> bool:
    text = value.strip()
    if not text:
        return True
    return bool(re.match(r"^(LOC\d+|A\d+N\d+_[a-z]+\d+|TRNA|TRN|rpl\d+|rps\d+)$", text, re.IGNORECASE))


def description_candidates(description: str) -> list[str]:
    text = clean_text(description)
    if not text:
        return []
    text = re.sub(r"[\[\]\(\),]", " ", text)
    parts = re.split(r"\s+", text)
    candidates: list[str] = []
    for part in parts:
        token = part.strip("\"' ;:")
        if len(token) < 3:
            continue
        if token.lower() in {"protein", "probable", "putative", "like", "family", "containing", "domain"}:
            continue
        if re.fullmatch(r"[A-Z0-9\-]{3,15}", token):
            candidates.append(token)
    return candidates


def choose_query_candidates(row: pd.Series) -> list[str]:
    candidates: list[str] = []
    gene = clean_text(row.get("gene", ""))
    description = clean_text(row.get("description", ""))
    gene_id = clean_text(row.get("gene_id", ""))

    for value in (gene, description, gene_id):
        value = clean_text(value)
        if not value:
            continue
        if value not in candidates and not looks_like_locus(value):
            candidates.append(value)

    for token in description_candidates(description):
        if token not in candidates:
            candidates.append(token)

    return candidates[:8]


def build_seed_query(row: pd.Series) -> str:
    candidates = choose_query_candidates(row)
    if candidates:
        return candidates[0]
    for field in ("gene", "description", "gene_id"):
        value = clean_text(row.get(field, ""))
        if value:
            return value
    return ""


def require_package(import_name: str, package_name: str):
    try:
        return __import__(import_name)
    except ImportError as exc:
        raise SystemExit(
            f"Missing dependency '{package_name}'. Install it before running this script."
        ) from exc


def map_with_mygene(frame: pd.DataFrame) -> pd.DataFrame:
    mygene = require_package("mygene", "mygene")
    mg = mygene.MyGeneInfo()

    mapped_rows: list[dict[str, object]] = []
    for _, row in frame.iterrows():
        row_dict = row.to_dict()
        candidates = choose_query_candidates(row)
        best_symbol = ""
        best_name = ""
        best_query = build_seed_query(row)
        best_score = math.nan
        best_raw_score = math.nan

        for candidate in candidates:
            try:
                results = mg.query(
                    candidate,
                    species="arabidopsis thaliana",
                    size=5,
                    fields="symbol,name,taxid,_score,alias",
                )
            except Exception:
                continue

            hits = [hit for hit in results.get("hits", []) if hit.get("symbol")]
            if not hits:
                continue

            hit = hits[0]
            score = hit.get("_score")
            if pd.notna(best_score) and pd.notna(score) and float(score) <= float(best_score):
                continue

            best_query = candidate
            best_symbol = str(hit.get("symbol", "")).strip()
            best_name = str(hit.get("name", "")).strip()
            best_score = score if score is not None else math.nan
            best_raw_score = score if score is not None else math.nan

        row_dict["ortholog_query"] = best_query
        row_dict["arabidopsis_symbol"] = best_symbol
        row_dict["arabidopsis_name"] = best_name
        row_dict["mapping_score"] = best_score
        row_dict["mapping_raw_score"] = best_raw_score
        mapped_rows.append(row_dict)

    return pd.DataFrame(mapped_rows)


def load_or_build_map(input_path: Path, output_path: Path, reuse_existing_maps: bool, skip_remap: bool) -> pd.DataFrame:
    if reuse_existing_maps and output_path.exists():
        return pd.read_csv(output_path)

    if skip_remap:
        if output_path.exists():
            return pd.read_csv(output_path)
        raise SystemExit(f"Cannot skip remap because {output_path} does not exist.")

    frame = pd.read_csv(input_path)
    mapped = map_with_mygene(frame)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    mapped.to_csv(output_path, index=False)
    return mapped


def unique_symbols(frame: pd.DataFrame) -> list[str]:
    if "arabidopsis_symbol" not in frame.columns:
        return []
    symbols = frame["arabidopsis_symbol"].fillna("").astype(str).str.strip()
    return sorted({symbol for symbol in symbols if symbol})


def run_gprofiler(symbols: Iterable[str]) -> pd.DataFrame:
    gprofiler = require_package("gprofiler", "gprofiler-official")
    gp = gprofiler.GProfiler(return_dataframe=True)
    symbol_list = list(symbols)
    if not symbol_list:
        return pd.DataFrame()
    result = gp.profile(
        organism="athaliana",
        query=symbol_list,
        sources=["GO:BP", "GO:MF", "GO:CC", "KEGG", "WP"],
    )
    if result is None:
        return pd.DataFrame()
    return pd.DataFrame(result)


def save_gprofiler_results(frame: pd.DataFrame, output_path: Path) -> pd.DataFrame:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(output_path, index=False)
    return frame


def plot_top_terms(frame: pd.DataFrame, output_path: Path, title: str, top_n: int) -> None:
    matplotlib = require_package("matplotlib", "matplotlib")
    plt = matplotlib.pyplot

    if frame.empty:
        fig, ax = plt.subplots(figsize=(8, 4))
        ax.text(0.5, 0.5, "No enrichment results", ha="center", va="center")
        ax.axis("off")
        fig.tight_layout()
        fig.savefig(output_path, dpi=200)
        plt.close(fig)
        return

    plot_frame = frame.sort_values("p_value").head(top_n).copy()
    plot_frame["minus_log10_p"] = -plot_frame["p_value"].clip(lower=1e-300).map(math.log10)
    plot_frame = plot_frame.iloc[::-1]

    fig_height = max(4, 0.45 * len(plot_frame) + 1.5)
    fig, ax = plt.subplots(figsize=(10, fig_height))
    ax.barh(plot_frame["name"], plot_frame["minus_log10_p"], color="#4C7A5A")
    ax.set_xlabel("-log10(p-value)")
    ax.set_title(title)
    ax.grid(axis="x", alpha=0.2)
    fig.tight_layout()
    fig.savefig(output_path, dpi=200)
    plt.close(fig)


def summary_dict(interaction_map: pd.DataFrame, flowering_map: pd.DataFrame) -> dict[str, int]:
    return {
        "significant_interaction_genes": int(len(interaction_map)),
        "mapped_interaction_genes": int((interaction_map["arabidopsis_symbol"].fillna("") != "").sum()),
        "mapped_interaction_symbols": int(len(unique_symbols(interaction_map))),
        "mapped_flowering_genes": int((flowering_map["arabidopsis_symbol"].fillna("") != "").sum()),
        "mapped_flowering_symbols": int(len(unique_symbols(flowering_map))),
    }


def markdown_table(frame: pd.DataFrame, top_n: int) -> str:
    if frame.empty:
        return "No enriched terms returned."

    cols = ["source", "native", "name", "p_value", "intersection_size", "term_size"]
    table = frame.loc[:, [col for col in cols if col in frame.columns]].sort_values("p_value").head(top_n).copy()
    if "p_value" in table.columns:
        table["p_value"] = table["p_value"].map(lambda value: f"{value:.6g}")
    return table.to_markdown(index=False)


def build_report(summary: dict[str, int], all_results: pd.DataFrame, flowering_results: pd.DataFrame, top_n: int) -> str:
    lines = [
        "# Ortholog-Based Enrichment for limma Interaction Genes",
        "",
        "This report maps hemp interaction genes to Arabidopsis ortholog symbols using local annotation text plus MyGene lookup, then runs g:Profiler enrichment in Arabidopsis space.",
        "The enrichment uses Arabidopsis genome-wide background, not a custom hemp transcriptome universe, so the results should be treated as pathway-level interpretation rather than exact odds-ratio calibration.",
        "",
        "## Mapping summary",
        "",
        f"- Significant limma interaction genes: {summary['significant_interaction_genes']}",
        f"- Interaction genes mapped to Arabidopsis-like orthologs: {summary['mapped_interaction_genes']}",
        f"- Unique Arabidopsis ortholog symbols from all interaction genes: {summary['mapped_interaction_symbols']}",
        f"- Flowering-focused interaction genes mapped: {summary['mapped_flowering_genes']}",
        f"- Unique Arabidopsis ortholog symbols from flowering-focused interaction genes: {summary['mapped_flowering_symbols']}",
        "",
        "## Top enriched terms from all interaction genes",
        "",
        markdown_table(all_results, top_n),
        "",
        "## Top enriched terms from flowering-focused interaction genes",
        "",
        markdown_table(flowering_results, top_n),
        "",
        "## Interpretation",
        "",
        "1. If plastid, thylakoid, and photosynthesis terms dominate the full interaction set, that indicates flowering stage differences are strongly coupled to cultivar-specific chloroplast and energy-state remodeling.",
        "2. If circadian rhythm, developmental timing, hormone signaling, or RNA-splicing terms appear in the flowering-focused subset, those are the strongest candidates for causal flowering regulation rather than downstream metabolic state changes.",
        "3. Use the flowering-focused mapped gene table first for breeding prioritization, and treat the broader chloroplast-enriched set as context for physiology and source-sink transition during floral induction.",
        "",
    ]
    return "\n".join(lines)


def main() -> None:
    args = parse_args()
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    interaction_map = load_or_build_map(
        INTERACTION_INPUT,
        INTERACTION_MAP,
        reuse_existing_maps=args.reuse_existing_maps,
        skip_remap=args.skip_remap,
    )
    flowering_map = load_or_build_map(
        FLOWERING_INPUT,
        FLOWERING_MAP,
        reuse_existing_maps=args.reuse_existing_maps,
        skip_remap=args.skip_remap,
    )

    all_results = run_gprofiler(unique_symbols(interaction_map))
    flowering_results = run_gprofiler(unique_symbols(flowering_map))

    save_gprofiler_results(all_results, ALL_GPROFILER)
    save_gprofiler_results(flowering_results, FLOWERING_GPROFILER)

    plot_top_terms(all_results, ALL_PLOT, "All interaction genes", args.top_n)
    plot_top_terms(flowering_results, FLOWERING_PLOT, "Flowering-focused interaction genes", args.top_n)

    summary = summary_dict(interaction_map, flowering_map)
    SUMMARY_JSON.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    REPORT_MD.write_text(build_report(summary, all_results, flowering_results, args.top_n), encoding="utf-8")


if __name__ == "__main__":
    main()