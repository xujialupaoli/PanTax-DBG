#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# 文件名: themis_scripts/tax_profile_strain_update.py

from pathlib import Path

def run(tax_profile, strain_abundance, output_tre):
    """
    根据计算出的菌株丰度 (strain_abundance)，更新原始 tax_profile.tre 中的丰度信息，
    并过滤掉丰度为 0 的分支，生成一个新的 tax_profile_strain.tre。
    """
    tax_file_path = Path(tax_profile)
    abundance_file_path = Path(strain_abundance)
    output_file_path = Path(output_tre)

    if not tax_file_path.exists():
        raise FileNotFoundError(f"Base tax profile not found: {tax_file_path}")
        
    if not abundance_file_path.exists():
        raise FileNotFoundError(f"Strain abundance file not found: {abundance_file_path}")
        

    print(f"[PanTax-DBG] Updating strain TRE...")
    print(f"  - Base Tree: {tax_file_path.name}")
    print(f"  - Strain Abundance: {abundance_file_path.name}")

    # 1. 读取菌株丰度表
    # 存储结构: {'taxid': 0.1, 'name': 0.1}
    strain_abundance_dict = {}
    # 存储原始字符串格式以保持精度: {'taxid': "0.1000..."}
    original_abundance_format = {}

    try:
        with open(abundance_file_path, 'r') as f:
            header = next(f, None) # 跳过表头
            for line in f:
                parts = line.strip().split('\t')
                # 兼容不同格式，通常是: TaxID, Name, Abundance 或 TaxID, Abundance
                if len(parts) >= 2:
                    # 假设第一列是ID，最后一列是丰度
                    s_id = parts[0].strip()
                    abund_str = parts[-1].strip()
                    try:
                        val = float(abund_str)
                        strain_abundance_dict[s_id] = val
                        original_abundance_format[s_id] = abund_str
                        # 如果有第二列是名称，也存一下
                        if len(parts) >= 3:
                            s_name = parts[1].strip()
                            strain_abundance_dict[s_name] = val
                            original_abundance_format[s_name] = abund_str
                    except ValueError:
                        continue
    except Exception as e:
        print(f"[PanTax-DBG] Error reading strain abundance: {e}")
        return

    # 2. 预读取分类树结构，建立 Strain -> Lineage 的映射
    # 我们需要知道哪些非strain行（如属、科）包含哪些有效的strain
    id_to_total_abundance = {} # 存储非strain层级的累计丰度
    
    try:
        with open(tax_file_path, 'r') as f_tax:
            lines = f_tax.readlines()

        # 第一遍遍历：计算所有上级节点的丰度
        for line in lines:
            parts = line.strip().split('\t')
            if len(parts) < 3: continue
            
            level = parts[0].strip().lower()
            current_id = parts[1].strip()
            # lineage通常在第三列，用|分隔
            lineage_ids = parts[2].strip().split('|')

            # 如果这一行是 strain，且我们在丰度表里有它
            if level == 'strain':
                # 尝试通过ID或名称获取丰度
                val = strain_abundance_dict.get(current_id, 0.0)
                
                # 如果这个strain有丰度，将其加到它路径上的所有父节点
                if val > 0:
                    for ancestor_id in lineage_ids:
                        id_to_total_abundance[ancestor_id] = id_to_total_abundance.get(ancestor_id, 0.0) + val
                    # 别忘了把自己也存进去
                    id_to_total_abundance[current_id] = val

        # 3. 第二遍遍历：写入新文件并过滤
        with open(output_file_path, 'w') as f_out:
            for line in lines:
                parts = line.strip().split('\t')
                if len(parts) < 3:
                    f_out.write(line)
                    continue

                current_id = parts[1].strip()
                
                # 获取该节点的计算丰度
                calculated_abundance = id_to_total_abundance.get(current_id, 0.0)
                
                # 过滤阈值 (1e-10)
                if calculated_abundance < 1e-10:
                    continue

                # 格式化丰度字符串
                # 如果是 Strain 层级且有原始字符串，优先用原始的
                if parts[0].strip().lower() == 'strain' and current_id in original_abundance_format:
                     final_abund_str = original_abundance_format[current_id]
                else:
                    # 避免科学计数法，保留10位小数并去除末尾0
                    final_abund_str = f"{calculated_abundance:.10f}".rstrip('0').rstrip('.')
                    if final_abund_str == '': final_abund_str = "0"

                # 替换原始行中的丰度列 (假设丰度在第9列/index 8，或者追加在最后)
                # PanTax-DBG/DBG-ganon 标准 TRE 格式：
                # rank, taxid, lineage, name, unique, shared, children, cum_abund, own_abund
                # 我们这里简化处理：保留前4列结构信息，追加/替换最后的新丰度
                
                # 如果原始行已经有足够的列，我们替换最后一列为新丰度(作为 relative abundance)
                # 或者按照您的需求，简单地追加到行尾，或者重构行
                # 这里采用最稳妥的方式：保留前8列结构，更新第9列(如果存在)或追加
                
                # 根据您提供的head数据：
                # unclassified - - unclassified 0 0 0 178061 35.61291
                # 这里的最后一列是相对丰度
                
                output_parts = parts[:]
                # 确保列表足够长
                while len(output_parts) < 9:
                    output_parts.append("0")
                
                # 更新最后一列 (Relative Abundance)
                output_parts[-1] = final_abund_str
                
                f_out.write("\t".join(output_parts) + "\n")
                
        print(f"[PanTax-DBG] Created strain-specific TRE: {output_file_path}")

    except Exception as e:
        print(f"[PanTax-DBG] Error writing strain TRE: {e}")
        import traceback
        traceback.print_exc()