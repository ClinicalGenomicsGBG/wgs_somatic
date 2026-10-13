def find_overlapping_gene(refseq_dict, variant_chrom, variant_pos):
    overlapping_list = []
    for gene in refseq_dict[variant_chrom]:
        for transcript in refseq_dict[variant_chrom][gene]:
            overlapping_dict = {}
            smallest = int(refseq_dict[variant_chrom][gene][transcript]["smallest"])
            largest = int(refseq_dict[variant_chrom][gene][transcript]["largest"])
            if variant_pos >= smallest and variant_pos <= largest:
                # variant is within transcript
                # within or between exons?
                for exon in refseq_dict[variant_chrom][gene][transcript]["exon"]:
                    start = int(
                        refseq_dict[variant_chrom][gene][transcript]["exon"][exon][
                            "start"
                        ]
                    )
                    stop = int(
                        refseq_dict[variant_chrom][gene][transcript]["exon"][exon][
                            "stop"
                        ]
                    )
                    if variant_pos >= start and variant_pos <= stop:
                        # variant is within exon
                        distance_start = abs(start - variant_pos)
                        distance_stop = abs(stop - variant_pos)
                        overlapping_dict["gene"] = gene
                        overlapping_dict["transcript"] = transcript
                        overlapping_dict["exon"] = f"exon{exon}"
                        overlapping_dict["minus_distance"] = distance_start
                        overlapping_dict["plus_distance"] = distance_stop
                        gene_info_string = f"{gene}:{transcript}:exon{exon} (-{distance_start}:+{distance_stop})"
                        overlapping_list.append(gene_info_string)
                        break
                    else:
                        if variant_pos < start:
                            # variant is between this exon and the previous one
                            previous_exon_stop = int(
                                refseq_dict[variant_chrom][gene][transcript]["exon"][
                                    previousexon
                                ]["stop"]
                            )
                            distance_upcoming_exon_start = abs(start - variant_pos)
                            distance_previous_exon_stop = abs(
                                previous_exon_stop - variant_pos
                            )
                            overlapping_dict["gene"] = gene
                            overlapping_dict["transcript"] = transcript
                            overlapping_dict["exon"] = f"exon{previousexon}-exon{exon}"
                            overlapping_dict["minus_distance"] = (
                                distance_previous_exon_stop
                            )
                            overlapping_dict["plus_distance"] = (
                                distance_upcoming_exon_start
                            )
                            gene_info_string = f"{gene}:{transcript}:exon{previousexon}(-{distance_previous_exon_stop})-exon{exon}(+{distance_upcoming_exon_start})"
                            overlapping_list.append(gene_info_string)
                            break
                        previousexon = exon
    return overlapping_list


def find_crossing_genes(refseq_dict, variant):
    variant_type = variant["INFO"]["SVTYPE"]
    if variant_type in ["CNV", "DEL", "DUP"]:
        gene_crossing = []
        v_chrom = variant["#CHROM"]
        start = int(variant["POS"])

        # Handle missing or invalid END field
        end_value = variant["INFO"].get("END", None)
        try:
            stop = int(end_value)
        except (ValueError, TypeError):
            # If END is missing or invalid, calculate it for deletions using SVLEN
            if variant_type == "DEL" and "SVLEN" in variant["INFO"]:
                svlen = int(variant["INFO"]["SVLEN"])
                stop = start + abs(svlen)
            else:
                print(
                    f"Warning: Missing or invalid END for variant {variant['ID']} at {v_chrom}:{start}"
                )
                stop = start  # Default to POS if END is missing or invalid

        try:
            for gene in refseq_dict[v_chrom]:
                smallest_list = []
                largest_list = []
                for transcript in refseq_dict[v_chrom][gene]:
                    g_start = int(refseq_dict[v_chrom][gene][transcript]["smallest"])
                    g_stop = int(refseq_dict[v_chrom][gene][transcript]["largest"])
                    smallest_list.append(g_start)
                    largest_list.append(g_stop)
                g_start = min(smallest_list)
                g_stop = max(largest_list)
                if start <= g_start:
                    if stop >= g_stop:
                        gene_crossing.append(
                            gene
                        )  # -------VS--GS---GE----VE INCLUDING GENE
                    elif stop >= g_start:
                        gene_crossing.append(
                            f"in:{gene}"
                        )  # ---VS---GS---VE--GE--- # INSIDE GENE
                    else:
                        pass  # ----VS---VE---GS---GE---- BEFORE GENE
                else:
                    if start <= g_stop:  # ------GS----VS--GE---
                        gene_crossing.append(
                            f"in:{gene}"
                        )  # ------GS--VS---VE--GE--- # INSIDE GENE
                    else:
                        pass  # ------GS---GE---VS---VE-- AFTER GENE
            gene_crossing_str = ",".join(gene_crossing)
            variant["gene_crossing"] = gene_crossing_str
        except Exception as e:
            print(f"Exception: {e}")
            variant["gene_crossing"] = None
    else:
        variant["gene_crossing"] = None
    return variant


"""
                # EXON PARSEING...Pain in the...brain
                if start <= g_start:
                    if stop >= g_stop:
                        variant["gene_crossing"].append(gene)  # -------VS--GS---GE----VE INCLUDING GENE
                    elif stop >= g_start:
                        # ------VS---GS---VE---GE INSIDE GENE
                        g_info = f"{gene}(exons:"
                        for exon in for refseq_dict[v_chrom][gene][transcript]["exon"]:
                            e_start = refseq_dict[v_chrom][gene][transcript]["exon"][exon]["start"]
                            e_stop = refseq_dict[v_chrom][gene][transcript]["exon"][exon]["stop"]
                            if stop >= e_stop:
                                g_info = f"{g_info}{exon}," # ------VS--GS--ES--EE--VE---GE INLCUDING EXON
                            elif stop >= e_start:
                                g_info = f"{g_info}{exon}:END" # ---VS--GS--ES--VE--EE--- INSIDE EXON
                            else:
                                pass # --------VS---GS--E1S--E1E--VE---E2S--E2E---GE-- OUTSIDE EXON
                        variant["gene_crossing"].append(g_info)
                    else:
                        pass
                        #----VS---VE---GS---GE---- BEFORE GENE
                else:
                    if start <= g_stop:
                        # --------GS---VS--GE
                        if stop <= g_stop:
                            # ----GS----VS---VE---GE
                            g_info = f"{gene}(exons:"
                            for exon in for refseq_dict[v_chrom][gene][transcript]["exon"]:
                                e_start = refseq_dict[v_chrom][gene][transcript]["exon"][exon]["start"]
                                e_stop = refseq_dict[v_chrom][gene][transcript]["exon"][exon]["stop"]
                                if start <= e_start: # ----GS---VS---ES----
                                    if stop >= e_stop:
                                        g_info = f"{g_info}{exon}," # -----GS---VS--ES--EE--VE---GE INCLUDING EXON
                                    elif stop >= e_start:
                                        g_info = f"{g_info}{exon}:END" # ---GS---VS--ES--VE--EE---GE--- INSIDE EXON
                                    else:
                                        # ---GS---VS--VE--EE--ES--GE--- OUTSIDE EXON
                                        pass
                                elif start >= e_start:
                                    # ----GS---ES---VS---
                                    if stop >=
                                    if start >= e_stop:
                                        

                                elif stop >= e_start:
                                    g_info = f"{g_info}{exon}:END"
                                else:
                                    pass # --------VS---GS--E1S--E1E--VE---E2S--E2E---GE
                            variant["gene_crossing"].append(g_info)
                        else:
                            # ----GS---VS--GE--VS--
                    else:
                        # ------GS---GE---VS---VE-- AFTER GENE
                    if stop <= g_stop:
                     start > g_start:
                    if start < g_stop:

        variant["gene_crossing"] = 
"""


def find_nearby_gene(refseq_dict, variant_chrom, variant_pos):
    upstream_distances = {}
    downstream_distances = {}
    for gene in refseq_dict[variant_chrom]:
        # find smallest upstream and downstream distances
        for transcript in refseq_dict[variant_chrom][gene]:
            smallest = int(refseq_dict[variant_chrom][gene][transcript]["smallest"])
            largest = int(refseq_dict[variant_chrom][gene][transcript]["largest"])
            if largest < variant_pos:
                # transcript is downstream of variant
                downstream_distances[f"{gene},{transcript}"] = abs(
                    variant_pos - largest
                )
            else:
                # transcript is upstream of variant
                upstream_distances[f"{gene},{transcript}"] = abs(variant_pos - smallest)
    try:
        min_upstream = min(upstream_distances, key=upstream_distances.get)
        up_gene, up_transcript = min_upstream.split(",")
        up_distance = upstream_distances[min_upstream]
    except:
        up_gene, up_transcript, up_distance = "N/A", "N/A", "N/A"

    try:
        min_downstream = min(downstream_distances, key=downstream_distances.get)
        down_gene, down_transcript = min_downstream.split(",")
        down_distance = downstream_distances[min_downstream]
    except:
        down_gene, down_transcript, down_distance = "N/A", "N/A", "N/A"

    nearby_dict = {}
    nearby_dict["down_gene"] = down_gene
    nearby_dict["down_transcript"] = down_transcript
    nearby_dict["down_distance"] = down_distance
    nearby_dict["up_gene"] = up_gene
    nearby_dict["up_transcript"] = up_transcript
    nearby_dict["up_distance"] = up_distance
    gene_info_string = [
        f"{down_gene}:{down_transcript} (-{down_distance}) | {up_gene}:{up_transcript} (+{up_distance})"
    ]

    return gene_info_string
