from cmath import nan
import enum
import pandas as pd
import re
import os


def manta_summary(mantaSV_vcf, mantaSV_summary, tumorname, genelist, normalname=''):
    
    df = pd.read_excel(str(mantaSV_vcf),engine='openpyxl')

    # Check for required columns
    required_columns = ['ID', 'DEL/DUP Genecrossings', 'GeneInfo 1', 'GeneInfo 2', 'FORMAT', tumorname]
    if normalname:
        required_columns.append(normalname)
    missing = [col for col in required_columns if col not in df.columns]
    if missing:
        print(f"Input Excel file is missing required columns: {missing}. Exiting manta_summary.")
        # Create an empty summary Excel with just headers
        with pd.ExcelWriter(str(mantaSV_summary)) as writer:
            df.to_excel(writer, sheet_name="Manta_raw", engine='xlsxwriter')
        return

    # drops rows with NA in columns ID. MATEID can be NA if not a translocation
    df = df.dropna(subset=['ID'])


    # this part of the script highlights genes from the gene list

    # open the genelist
    genelist = open(genelist, "r")

    genelist = genelist.readlines()

    # make a list with genes from the genelist
    gene_list = []
    for gene in genelist:
        gene_list.append(gene.rstrip("\n"))
    # remove duplicates 
    gene_list = list(set(gene_list))

    # pattern to search for with an "or" between each gene
    pattern = '|'.join(gene_list)
    # search for patterns of the genes in the genelist in certain column
    if df['DEL/DUP Genecrossings'].isnull().all():
        column_patterns_genecrossings = False
    else:
        column_patterns_genecrossings = df.loc[df['DEL/DUP Genecrossings'].str.contains(pattern, na=False)]['DEL/DUP Genecrossings']
    column_patterns_geneinfo1 = df.loc[df['GeneInfo 1'].str.contains(pattern, na=False)]['GeneInfo 1']
    column_patterns_geneinfo2 = df.loc[df['GeneInfo 2'].str.contains(pattern, na=False)]['GeneInfo 2']
    

    # create new column Genelist to contain genes from genelist that exist in certain columns
    df['Genelist'] = ""

    # function to find only whole words, otherwise it would find gene AR in XXARXX
    def find_only_whole_word(search_string, input_string):
        # Create a raw string with word boundaries from the user's input_string
        raw_search_string = r"\b" + search_string + r"\b"

        match_output = re.search(raw_search_string, input_string)

        no_match_was_found = ( match_output is None )
        if no_match_was_found:
            return False
        else:
            return True
    
    # function appending genes found to Genelist column
    def append_genes(col_patterns, col):
        for gene in gene_list:
            for row_of_genes in col_patterns:
                if gene in row_of_genes:
                    
                    if find_only_whole_word(gene, row_of_genes) == True:
                        
                        if len(df[col==row_of_genes].index.values) > 1:
                            
                            for rowno in df[col==row_of_genes].index.values:
                                
                                index_value=int(rowno)
                                df.at[index_value, 'Genelist'] = df.at[index_value, 'Genelist'] + gene + ' '

                        else:
                            index_value = int(df[col==row_of_genes].index.values)
                            df.at[index_value, 'Genelist'] = df.at[index_value, 'Genelist'] + gene + ' '
        return df

    if not column_patterns_genecrossings is False:
        append_genes(column_patterns_genecrossings, df['DEL/DUP Genecrossings'])
    append_genes(column_patterns_geneinfo1, df['GeneInfo 1'])
    append_genes(column_patterns_geneinfo2, df['GeneInfo 2'])

    
    # remove duplicates in Genelist column
    for ind in df.index: 
        df['Genelist'][ind] = ' '.join(set(df['Genelist'][ind].split()))


    # this part of the script extracts allele frequencies from PR/SR

    # extract variants only supported by both PR and SR
    df = df.loc[df['FORMAT'] == 'PR:SR']

    if normalname:  
        # add columns for PR/SR for normal sample
        df[normalname + ':PR'] = ""
        df[normalname + ':PR-alt'] = ""
        df[normalname + ':SR'] = ""
        df[normalname + ':SR-alt'] = ""
        df['TOTAL alt (N)'] = ""
        df['TOTAL VAF (N)'] = ""

        # for row in column that contains PR/SR info for normal, add info to new columns
        for row in df[normalname]:
            PR = int(re.split('[,:]', row)[0])
            PR_alt = int(re.split('[,:]', row)[1])
            SR = int(re.split('[,:]', row)[2])
            SR_alt = int(re.split('[,:]', row)[3])
        
            if len(df.loc[df[normalname] == row].index.values) > 1:
                for rown in df.loc[df[normalname] == row].index.values:
                    ind_val=int(rown)
                    df.at[ind_val, normalname + ':PR'] = PR
                    df.at[ind_val, normalname + ':PR-alt'] = PR_alt
                    df.at[ind_val, normalname + ':SR'] = SR
                    df.at[ind_val, normalname + ':SR-alt'] =  SR_alt
                    df.at[ind_val, 'TOTAL alt (N)'] = PR_alt + SR_alt
                    if PR + PR_alt + SR + SR_alt == 0:
                        df.at[ind_val, 'TOTAL VAF (N)'] = ''
                    else:
                        df.at[ind_val, 'TOTAL VAF (N)'] = str(int(round(float(PR_alt + SR_alt) / (PR + PR_alt + SR + SR_alt) *100))) + '%'        
            else:
                row_index = int(df.loc[df[normalname] == row].index.values)
                df.at[row_index, normalname + ':PR'] = PR
                df.at[row_index, normalname + ':PR-alt'] = PR_alt
                df.at[row_index, normalname + ':SR'] = SR
                df.at[row_index, normalname + ':SR-alt'] =  SR_alt
                df.at[row_index, 'TOTAL alt (N)'] = PR_alt + SR_alt

                if PR + PR_alt + SR + SR_alt == 0:
                    df.at[row_index, 'TOTAL VAF (N)'] = ''
                else:
                    df.at[row_index, 'TOTAL VAF (N)'] = str(int(round(float(PR_alt + SR_alt) / (PR + PR_alt + SR + SR_alt) *100))) + '%'
    

    # add columns for PR/SR for tumor sample
    df[tumorname + ':PR'] = ""
    df[tumorname + ':PR-alt'] = ""
    df[tumorname + ':SR'] = ""
    df[tumorname + ':SR-alt'] = ""
    df['TOTAL alt (T)'] = ""
    df['TOTAL VAF (T)'] = ""

    # for row in column that contains PR/SR info for tumor, add info to new columns
    for row in df[tumorname]:
        
        PR = int(re.split('[,:]', row)[0])
        PR_alt = int(re.split('[,:]', row)[1])
        SR = int(re.split('[,:]', row)[2])
        SR_alt = int(re.split('[,:]', row)[3])
        
        if len(df.loc[df[tumorname] == row].index.values) > 1:
            for rown in df.loc[df[tumorname] == row].index.values:
                ind_val=int(rown)
                df.at[ind_val, tumorname + ':PR'] = PR
                df.at[ind_val, tumorname + ':PR-alt'] = PR_alt
                df.at[ind_val, tumorname + ':SR'] = SR
                df.at[ind_val, tumorname + ':SR-alt'] =  SR_alt
                df.at[ind_val, 'TOTAL alt (T)'] = PR_alt + SR_alt
                if PR + PR_alt + SR + SR_alt == 0:
                    df.at[ind_val, 'TOTAL VAF (T)'] = ''
                else:
                    df.at[ind_val, 'TOTAL VAF (T)'] = str(int(round(float(PR_alt + SR_alt) / (PR + PR_alt + SR + SR_alt) *100))) + '%'        
        else:
            row_index = int(df.loc[df[tumorname] == row].index.values)
            df.at[row_index, tumorname + ':PR'] = PR
            df.at[row_index, tumorname + ':PR-alt'] = PR_alt
            df.at[row_index, tumorname + ':SR'] = SR
            df.at[row_index, tumorname + ':SR-alt'] =  SR_alt
            df.at[row_index, 'TOTAL alt (T)'] = PR_alt + SR_alt
            if PR + PR_alt + SR + SR_alt == 0:
                df.at[row_index, 'TOTAL VAF (T)'] = ''
            else:
                df.at[row_index, 'TOTAL VAF (T)'] = str(int(round(float(PR_alt + SR_alt) / (PR + PR_alt + SR + SR_alt) *100))) + '%'

    # Second df with a selection of columns

    if normalname:
        df2 = df[["Varianttype", "Breakpoint 1", "GeneInfo 1", "Breakpoint 2", "GeneInfo 2", "ALT", "FORMAT", normalname, "TOTAL VAF (N)", tumorname, "TOTAL VAF (T)", "DEL/DUP Genecrossings", "Genelist"]]
    else:
        df2 = df[["Varianttype", "Breakpoint 1", "GeneInfo 1", "Breakpoint 2", "GeneInfo 2", "ALT", "FORMAT", tumorname, "TOTAL VAF (T)", "DEL/DUP Genecrossings", "Genelist"]]

    # If more than 30 genes in DEL/DUP Genecrossings - write number of genes instead of names of genes 
    for deldup in df2['DEL/DUP Genecrossings']:
        deldup = str(deldup)
        in_genes = []
        for gene in deldup.split(','):
            if 'in:' in gene:
                in_genes.append(gene)
        if len(deldup.split(',')) > 30:
            no_of_genes = len(deldup.split(','))
            new_deldup = f'No of genes: {no_of_genes}. \n {", ".join(in_genes)}'
            df2['DEL/DUP Genecrossings'] = df2['DEL/DUP Genecrossings'].replace(deldup, new_deldup)


    with pd.ExcelWriter(str(mantaSV_summary)) as writer:
        # Write both dfs to same excel file but different sheets
        df.to_excel(writer, sheet_name="Manta_raw", engine='xlsxwriter')
        df2.to_excel(writer, sheet_name="Manta_report", engine='xlsxwriter')

        # Modify Manta_report sheet
        workbook = writer.book
        worksheet = writer.sheets["Manta_report"]
        wrap_format = workbook.add_format({'text_wrap': True, 'font_name': 'Cambria (Headings)', 'font_size': 14})

        # Adjust cell width of column
        for column in df2:
            col_idx = df2.columns.get_loc(column)
            worksheet.set_column(col_idx+1, col_idx+1, len(column)+10)

        # Adjust row height
        row_idx = 1
        for index, row in df2.iterrows():
            column_idx = 0
            for cell in row:
                column_idx +=1
                try:
                    worksheet.write(row_idx, column_idx, cell, wrap_format)
                except:
                    continue
            row_idx += 1
    

        # Adjust font, color, etc
        header_format = workbook.add_format()
        header_format.set_font_name('Cambria (Headings)')
        header_format.set_font_size(16)
        header_format.set_bg_color('#99CCFF')
        header_format.set_bold()

        # Write adjustments to sheet
        for col_num, value in enumerate(df2.columns.values):
            worksheet.write(0, col_num+1, value, header_format)


        writer.close()
