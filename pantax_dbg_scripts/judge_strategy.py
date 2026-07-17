# -*- coding: utf-8 -*-

def get_top_two_abundance_from_list(val_list):
    if not val_list:
        return 0.0, 0.0
    sorted_vals = sorted(val_list, reverse=True)
    top1 = sorted_vals[0]
    top2 = sorted_vals[1] if len(sorted_vals) > 1 else 0.0
    return top1, top2

def calculate_decision(top1_ganon, top2_ganon, top1_ggcat, top2_ggcat):
    """Select an abundance-allocation strategy from the two leading estimates."""

    # Relative separation between the two most abundant species.
    rel_diff_ganon = (top1_ganon - top2_ganon) / top1_ganon if top1_ganon > 1e-9 else 0.0
    rel_diff_ggcat = (top1_ggcat - top2_ggcat) / top1_ggcat if top1_ggcat > 1e-9 else 0.0
    
    i = rel_diff_ggcat - rel_diff_ganon

    is_low_abundance = (top1_ggcat < 0.12) and (top1_ganon < 0.12)

    # Guard against large cross-strategy discrepancies in low-abundance profiles.
    if is_low_abundance:
        safe_top1_ganon = top1_ganon if top1_ganon > 1e-9 else 1e-9
        safe_top1_ggcat = top1_ggcat if top1_ggcat > 1e-9 else 1e-9

        if (top1_ggcat / safe_top1_ganon) > 2.0:
            return "GANON"

        if (top1_ganon / safe_top1_ggcat) > 2.0:
            return "GGCAT"

    # Select the strategy with the more suitable leading-rank contrast.
    if i < 0:
        strategy = "GANON"

        if is_low_abundance and i < -0.08:
            strategy = "GGCAT"

        return strategy

    else:
        strategy = "GGCAT"

        if is_low_abundance and i > 0.08:
            strategy = "GANON"

        return strategy
