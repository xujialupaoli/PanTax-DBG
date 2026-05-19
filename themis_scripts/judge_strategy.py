# -*- coding: utf-8 -*-

def get_top_two_abundance_from_list(val_list):
    """从丰度列表中获取前二名"""
    if not val_list:
        return 0.0, 0.0
    sorted_vals = sorted(val_list, reverse=True)
    top1 = sorted_vals[0]
    top2 = sorted_vals[1] if len(sorted_vals) > 1 else 0.0
    return top1, top2

def calculate_decision(top1_ganon, top2_ganon, top1_ggcat, top2_ggcat):
    """
    根据相对差异 (i) 和 Top1 比值 (Ratio) 决定权重策略。
    
    逻辑:
    1. 计算 i = RD_GGCAT - RD_GANON
    2. 判定是否为低丰度场景 (<0.12)
    3. 低丰度下，如果 Top1 比值 > 2，强制选小的那个。
    4. 否则根据 i 的正负判定，并剔除低丰度下的反常 i 值 (<-0.08 或 >0.08)。
    """
    
    # 1. 计算相对差异 (Relative Difference)
    # 避免除以 0
    rel_diff_ganon = (top1_ganon - top2_ganon) / top1_ganon if top1_ganon > 1e-9 else 0.0
    rel_diff_ggcat = (top1_ggcat - top2_ggcat) / top1_ggcat if top1_ggcat > 1e-9 else 0.0
    
    # 2. 计算指标 i
    i = rel_diff_ggcat - rel_diff_ganon
    
    # 3. 判定低丰度 (相对均衡)
    is_low_abundance = (top1_ggcat < 0.12) and (top1_ganon < 0.12)
    
    # =======================================================
    # 优先级 1: 低丰度下的 Top1 倍数悬殊检测 (Ratio Override)
    # =======================================================
    if is_low_abundance:
        # 避免除以 0
        safe_top1_ganon = top1_ganon if top1_ganon > 1e-9 else 1e-9
        safe_top1_ggcat = top1_ggcat if top1_ggcat > 1e-9 else 1e-9

        # 如果 GGCAT 的 Top1 是 GANON 的 2 倍以上 -> 认为 GGCAT 冒进 -> 用 GANON
        if (top1_ggcat / safe_top1_ganon) > 2.0:
            return "GANON"

        # 如果 GANON 的 Top1 是 GGCAT 的 2 倍以上 -> 认为 GANON 冒进 -> 用 GGCAT
        if (top1_ganon / safe_top1_ggcat) > 2.0:
            return "GGCAT"

    # =======================================================
    # 优先级 2: 基于 i 值的常规与反常检测
    # =======================================================
    
    # === 分支 A: i < 0 (默认 GANON) ===
    if i < 0:
        # 默认
        strategy = "GANON"
        
        # 反常检测: 低丰度下，如果 i 非常小 (<-0.08)，说明 Ganon 结果可能没那么好
        if is_low_abundance and i < -0.08:
            strategy = "GGCAT"
            
        return strategy

    # === 分支 B: i >= 0 (默认 GGCAT) ===
    else:
        # 默认
        strategy = "GGCAT"
        
        # 反常检测: 低丰度下，如果 i 非常大 (>0.08)，说明 GGCAT 结果可能没那么好
        if is_low_abundance and i > 0.08:
            strategy = "GANON"
            
        return strategy