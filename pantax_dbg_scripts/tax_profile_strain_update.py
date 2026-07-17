#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from pathlib import Path

def run(tax_profile, strain_abundance, output_tre):
    """Replace tree abundances with sums from the final strain profile."""
    tax_file_path = Path(tax_profile)
    abundance_file_path = Path(strain_abundance)
    output_file_path = Path(output_tre)

    if not tax_file_path.exists():
        raise FileNotFoundError(f"Base tax profile not found: {tax_file_path}")
        
    if not abundance_file_path.exists():
        raise FileNotFoundError(f"Strain abundance file not found: {abundance_file_path}")
        

    print("[PanTax-DBG] Updating strain TRE...")
    print(f"  - Base Tree: {tax_file_path.name}")
    print(f"  - Strain Abundance: {abundance_file_path.name}")

    # Index strain abundances by both taxonomic identifier and strain name.
    strain_abundance_dict = {}
    original_abundance_format = {}

    try:
        with open(abundance_file_path, 'r') as f:
            next(f, None)
            for line in f:
                parts = line.strip().split('\t')
                # Accept either taxid-name-abundance or taxid-abundance rows.
                if len(parts) >= 2:
                    s_id = parts[0].strip()
                    abund_str = parts[-1].strip()
                    try:
                        val = float(abund_str)
                        strain_abundance_dict[s_id] = val
                        original_abundance_format[s_id] = abund_str
                        if len(parts) >= 3:
                            s_name = parts[1].strip()
                            strain_abundance_dict[s_name] = val
                            original_abundance_format[s_name] = abund_str
                    except ValueError:
                        continue
    except Exception as e:
        print(f"[PanTax-DBG] Error reading strain abundance: {e}")
        return

    # Accumulate each strain abundance along its taxonomic lineage.
    id_to_total_abundance = {}
    
    try:
        with open(tax_file_path, 'r') as f_tax:
            lines = f_tax.readlines()

        for line in lines:
            parts = line.strip().split('\t')
            if len(parts) < 3: continue
            
            level = parts[0].strip().lower()
            current_id = parts[1].strip()
            lineage_ids = parts[2].strip().split('|')

            if level == 'strain':
                val = strain_abundance_dict.get(current_id, 0.0)

                if val > 0:
                    for ancestor_id in lineage_ids:
                        id_to_total_abundance[ancestor_id] = id_to_total_abundance.get(ancestor_id, 0.0) + val
                    id_to_total_abundance[current_id] = val

        # Retain positive-abundance branches and update the final tree column.
        with open(output_file_path, 'w') as f_out:
            for line in lines:
                parts = line.strip().split('\t')
                if len(parts) < 3:
                    f_out.write(line)
                    continue

                current_id = parts[1].strip()
                
                calculated_abundance = id_to_total_abundance.get(current_id, 0.0)

                if calculated_abundance < 1e-10:
                    continue

                if parts[0].strip().lower() == 'strain' and current_id in original_abundance_format:
                     final_abund_str = original_abundance_format[current_id]
                else:
                    final_abund_str = f"{calculated_abundance:.10f}".rstrip('0').rstrip('.')
                    if final_abund_str == '': final_abund_str = "0"

                
                
                output_parts = parts[:]
                while len(output_parts) < 9:
                    output_parts.append("0")
                
                output_parts[-1] = final_abund_str
                
                f_out.write("\t".join(output_parts) + "\n")
                
        print(f"[PanTax-DBG] Created strain-specific TRE: {output_file_path}")

    except Exception as e:
        print(f"[PanTax-DBG] Error writing strain TRE: {e}")
        import traceback
        traceback.print_exc()
