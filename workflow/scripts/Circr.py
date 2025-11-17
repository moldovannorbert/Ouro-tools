import pandas as pd
pd.options.mode.chained_assignment = None
import numpy as np
import sys,re,os,time,argparse,pybedtools
import subprocess
import logging
from functools import partial
from multiprocessing import Pool
from itertools import islice
import tempfile
import shutil
import concurrent.futures
from concurrent.futures import ThreadPoolExecutor
import traceback

#@#@#@# DICTIONARY SECTION #@#@#@#

match = {'CG': '|',
		 'GC': '|',
		 'AU': '|',
		 'UA': '|',
		 'UG': ':',
		 'GU': ':'}

seeds_a = {' |||||| ': '7A1',
	  	   ':|||||| ': '7mer-A1',
		   '||||||| ': '8mer',
		   '||| ||| ': 'SM',
		   '|||:||| ': 'SM',
		   '||||||  ': 'Off-6mer'}

seeds_b = {' |||||| ': '6mer',
		   ':|||||| ': '6mer',
		   '||||||| ': '7mer-m8'}

DataTree = {'human': {
				'file': 'support_files/miRNA/hsa_mature.fa',
				'acronym': 'hsa',
				'hg19':{
						'AGO': 'support_files/human/hg19/hg19.AGO.bed',
						'INT': 'support_files/human/hg19/hg19.INT.bed',
						'GTF': 'support_files/human/hg19/hg19.ensGene.gtf',
		 				'TAXA': 9606,
		 				'Genome': 'support_files/human/hg19/hg19.fa',
		 				'rRNA_Coords': 'support_files/human/hg19/hg19.rRNA.bed'
						},
		 		'hg38':{
						'AGO': 'support_files/human/hg38/hg38.AGO.bed',
						'INT': 'support_files/human/hg38/hg38.INT.bed',
 						'GTF': 'support_files/human/hg38/hg38.ensGene.gtf',
 		 				'TAXA': 9606,
 		 				'Genome': 'support_files/human/hg38/hg38.fa',
 		 				'rRNA_Coords': 'support_files/human/hg38/hg38.rRNA.bed'
		 				}
		  		},
			'mouse': {
				'file' : 'support_files/miRNA/mmu_mature.fa',
				'acronym' : 'mmu',
				'mm9':{
						'AGO': 'support_files/mouse/mm9/mm9.AGO.bed',
						'INT': 'support_files/mouse/mm9/mm9.INT.bed',
						'GTF': 'support_files/mouse/mm9/mm9.ensGene.gtf',
		 				'TAXA': 10090,
		 				'Genome': 'support_files/mouse/mm9/mm9.fa',
		 				'rRNA_Coords': 'support_files/mouse/mm9/mm9.rRNA.bed'
					  },
		 		'mm10':{
						'AGO': 'support_files/mouse/mm10/mm10.AGO.bed',
						'INT': 'support_files/mouse/mm10/mm10.INT.bed',
 						'GTF': 'support_files/mouse/mm10/mm10.ensGene.gtf',
 		 				'TAXA': 10090,
 		 				'Genome': 'support_files/mouse/mm10/mm10.fa',
 		 				'rRNA_Coords': 'support_files/mouse/mm10/mm10.rRNA.bed'
		 			   }
			     },
			'worm': {
				'file': 'support_files/miRNA/cel_mature.fa',
				'acronym': 'cel',
				'ce10':{
						'AGO': 'support_files/worm/ce10/ce10.AGO.bed',
						'INT': 'support_files/worm/ce10/ce10.INT.bed',
						'GTF': 'support_files/worm/ce10/ce10.ensGene.gtf',
		 				'TAXA': 6239,
		 				'Genome': 'support_files/worm/ce10/ce10.fa',
		 				'rRNA_Coords': 'support_files/worm/ce10/ce10.rRNA.bed'
					   },
			    'ce11':{
						'AGO': 'support_files/worm/ce11/ce11.AGO.bed',
						'INT': 'support_files/worm/ce11/ce11.INT.bed',
 						'GTF': 'support_files/worm/ce11/ce11.ensGene.gtf',
 		 				'TAXA': 6239,
 		 				'Genome': 'support_files/worm/ce11/ce11.fa',
 		 				'rRNA_Coords': 'support_files/worm/ce11/ce11.rRNA.bed'
		 			   }
				},
			'fruitfly': {
					'file': 'support_files/miRNA/dme_mature.fa',
					'acronym': 'dme',
					'dm3':{
							'AGO': 'support_files/fruitfly/dm3/dm3.AGO.bed',
							'INT': 'support_files/fruitfly/dcm3/dm3.INT.bed',
							'GTF': 'support_files/fruitfly/dm3/dm3.ensGene.gtf',
		 					'TAXA': 7227,
		 					'Genome': 'support_files/fruitfly/dm3/dm3.fa',
		 					'rRNA_Coords': 'support_files/fruitfly/dm3/dm3.rRNA.bed'
						  },
		    		'dm6':{
							'AGO': 'support_files/fruitfly/dm6/dm6.AGO.bed',
							'INT': 'support_files/fruitfly/dm6/dm6.INT.bed',
 							'GTF': 'support_files/fruitfly/dm6/dm6.ensGene.gtf',
 		 					'TAXA': 7227,
 		 					'Genome': 'support_files/fruitfly/dm6/dm6.fa',
 		 					'rRNA_Coords': 'support_files/fruitfly/dm6/dm6.rRNA.bed'
		 				  }
		   	  		}
		}

organisms = {'human':['hg19', 'hg38'],
			 'fruitfly':['dm3', 'dm6'],
			 'worm':['ce10', 'ce11'],
			 'mouse': ['mm9', 'mm10']}

#@#@#@# DEFINITION SECTION #@#@#@#

def getCircrPath():
	c_path = os.path.realpath(__file__).replace("Circr.py","")
	return c_path

def hybrid_parser(lines, DATA, key):
	df = pd.DataFrame(columns=range(len(lines[9].strip()[10:-3])))
	df.loc[len(df)] = [x for x in lines[9].strip()[10:-3]]
	df.loc[len(df)] = [x for x in lines[10].strip('\n')[10:-3]]
	df.loc[len(df)] = [x for x in lines[11].strip('\n')[10:-3]]
	df.loc[len(df)] = [x for x in lines[12].strip()[10:-3]]
	m = []
	target = df.head(2).max().to_list()
	mirna = df.tail(2).max().to_list()
	tuples = tuple(zip(target, mirna))
	for i in tuples:
		if ''.join(i) in match.keys():
			m.append(match[''.join(i)])
		else:
			m.append(' ')
	ll = m.index('|')
	rl = -(m[::-1].index('|')+1)
	try:
		if (-(m[::-1].index(':'))+1) > rl:
			rl = -(m[::-1].index(':')+1)
	except ValueError:
		pass
	upper = ''.join([''.join(target[:ll]).lower(), ''.join(target[ll:rl]).replace(' ','-'), ''.join(target[rl:]).lower()])
	lower = ''.join([''.join(mirna[:ll]).lower(), ''.join(mirna[ll:rl]).replace(' ','-'), ''.join(mirna[rl:]).lower()])
	if (upper.startswith('a')) and (''.join(m[-8:]) in seeds_a):
		DATA[key].append(seeds_a[''.join(m[-8:])])
	elif (not upper.startswith('a')) and (''.join(m[-8:]) in seeds_b):
		DATA[key].append(seeds_b[''.join(m[-8:])])
	else:
		DATA[key].append('No Seed')

def run_process(obj):
    try:
        result = subprocess.run(obj, shell=True, check=True, capture_output=True, text=True)
        if result.stdout:
            logging.debug(f"Process output: {result.stdout.strip()}")
        return True
    except subprocess.CalledProcessError as e:
        logging.error(f"Process failed with exit code {e.returncode}: {e.stderr.strip()}")
        return False
    except Exception as ex:
        logging.error(f"Error running process: {str(ex)}")
        return False

def examine(fileIn, DB, acronym, keep_tmp=False):
    # Get the directory and base filename
    file_dir = os.path.dirname(fileIn)
    base_name = os.path.basename(fileIn)
    
    # Construct the input file path
    mir_tmp_file = os.path.join(file_dir, 'mir_TMP_' + base_name)
    
    # Skip if the file doesn't exist
    if not os.path.exists(mir_tmp_file):
        logging.warning(f"miRanda output file {mir_tmp_file} not found. Skipping.")
        return DB
    
    DATA = {}
    key = ''
    processed_pairs = set()  # Track processed miRNA-circRNA pairs to avoid duplicates
    new_rows = []  # Collect all rows before adding to DB
    
    with open(mir_tmp_file, 'r') as file:
        for line in file:
            if "Performing" in line.strip().split(' '):
                key = ''.join(line.strip().split(' ')[2:])
                if key not in DATA.keys():
                    DATA[key] = []
            elif "Query:" in line.strip().split(' '):
                DATA[key].append(' '.join(line.strip().split(' ')[4:]))
            elif "|" in line.strip():
                DATA[key].append(line.strip('\n')[16:])
            elif "Ref:" in line.strip().split(' '):
                DATA[key].append(' '.join(line.strip().split(' ')[6:]))
            elif line.strip().split(' ')[0].startswith('>{0}'.format(acronym)):
                DATA[key].append(line.strip().split('\t')[2])
                DATA[key].append(line.strip().split('\t')[3])
                DATA[key].extend(line.strip().split('\t')[5].split(' '))
    
    # Filter out empty values
    DATA = {key: value for key, value in DATA.items() if value}
    
    # Process all data at once
    for key in DATA.keys():
        mirna_name = key.split('vs')[0]
        circ_name = key.split('vs')[1]
        
        # Check if this pair has been processed already
        mirna_circ_pair = (mirna_name, circ_name)
        if mirna_circ_pair not in processed_pairs:
            processed_pairs.add(mirna_circ_pair)
            
            if (DATA[key][2][-4] == 'a') and (DATA[key][1][-8:] in seeds_a):
                DATA[key].append(seeds_a[DATA[key][1][-8:]])
            elif (DATA[key][2][-4] != 'a') and (DATA[key][1][-8:] in seeds_b):
                DATA[key].append(seeds_b[DATA[key][1][-8:]])
            else:
                DATA[key].append('No Seed')
            
            # Extract chromosome and strand from Circ Name if available
            chrom = ''
            strand = ''
            try:
                # Format is typically chr:start-end or chr:start-end:strand
                circ_parts = circ_name.split(':')
                if len(circ_parts) >= 1:
                    chrom = circ_parts[0]
                
                # Try to extract strand from the end
                if len(circ_parts) > 2 and circ_parts[-1] in ['+', '-']:
                    strand = circ_parts[-1]
                else:
                    # Try to extract from a pattern like chr1:12345-67890(+)
                    strand_match = re.search(r'[()]([+-])[()]', circ_name)
                    if strand_match:
                        strand = strand_match.group(1)
                    else:
                        strand = '+'  # Default to + if not found
            except Exception as e:
                logging.warning(f"Error extracting chromosome/strand from {circ_name}: {e}")
            new_rows.append({
                'Chrom': chrom,
                'miRNA Name': mirna_name,
                'Circ Name': circ_name,
                'Strand': strand,
                'Seed Category': DATA[key][-1],
                'Start': DATA[key][-3],
                'End': DATA[key][-2],
                'ID': 'INT_M_'+str(len(DB) + len(new_rows))
            })
    
    # Add new rows to DB
    if new_rows:
        new_df = pd.DataFrame(new_rows)
        DB = pd.concat([DB, new_df], ignore_index=True)
    
    if not keep_tmp:
        os.remove(mir_tmp_file)

    return DB

def hybrid(fileIn, keep_tmp=False):
	# Get the directory and base filename
	file_dir = os.path.dirname(fileIn)
	base_name = os.path.basename(fileIn)
	
	# Construct the input and output file paths
	hyb_tmp_file = os.path.join(file_dir, 'hyb_TMP_' + base_name)
	hyb_parsed_file = os.path.join(file_dir, 'hyb_Parsed_' + base_name)
	
	# Skip processing if the parsed file already exists and keep_tmp is True
	if keep_tmp and os.path.exists(hyb_parsed_file):
		logging.debug(f"Skipping RNAhybrid parsing for {base_name} as parsed file already exists")
		return
	
	DF = pd.DataFrame(columns = ['miRNA Name','Circ Name', 'Seed Category', 'Start', 'End'])
	DATA = {}
	key = ''
	counter = 15
	processed_pairs = set()  # Track processed miRNA-circRNA pairs to avoid duplicates
	
	with open(hyb_tmp_file, 'r') as file:
		while True:
			try:
				lines = list(islice(file, counter))
				if lines[0].startswith('target too long'):
					del lines[0]
					counter = 16
					key = lines[0].strip().split(': ')[1] + 'vs' + lines[2].strip().split(': ')[1]
					# Check if this pair has been processed already
					mirna_circ_pair = (key.split('vs')[1], key.split('vs')[0])
					if mirna_circ_pair not in processed_pairs:
						if key not in DATA.keys():
							DATA[key] = [int(lines[8].strip().split(' ')[-1]),int(lines[3].strip().split(': ')[1])+int(lines[8].strip().split(' ')[-1])]
							hybrid_parser(lines, DATA, key)
							DF.loc[len(DF)] = [key.split('vs')[1], key.split('vs')[0], DATA[key][-1], DATA[key][0], DATA[key][1]]
							processed_pairs.add(mirna_circ_pair)
				elif lines[0].startswith('\n'):
					del lines[0]
					counter = 15
					key = lines[0].strip().split(': ')[1] + 'vs' + lines[2].strip().split(': ')[1]
					# Check if this pair has been processed already
					mirna_circ_pair = (key.split('vs')[1], key.split('vs')[0])
					if mirna_circ_pair not in processed_pairs:
						if key not in DATA.keys():
							DATA[key] = [int(lines[8].strip().split(' ')[-1]),int(lines[3].strip().split(': ')[1])+int(lines[8].strip().split(' ')[-1])]
							hybrid_parser(lines, DATA, key)
							DF.loc[len(DF)] = [key.split('vs')[1], key.split('vs')[0], DATA[key][-1], DATA[key][0], DATA[key][1]]
							processed_pairs.add(mirna_circ_pair)
				elif lines[0].startswith('target:'):
					counter = 15
					key = lines[0].strip().split(': ')[1] + 'vs' + lines[2].strip().split(': ')[1]
					# Check if this pair has been processed already
					mirna_circ_pair = (key.split('vs')[1], key.split('vs')[0])
					if mirna_circ_pair not in processed_pairs:
						if key not in DATA.keys():
							DATA[key] = [int(lines[8].strip().split(' ')[-1]),int(lines[3].strip().split(': ')[1])+int(lines[8].strip().split(' ')[-1])]
							hybrid_parser(lines, DATA, key)
							DF.loc[len(DF)] = [key.split('vs')[1], key.split('vs')[0], DATA[key][-1], DATA[key][0], DATA[key][1]]
							processed_pairs.add(mirna_circ_pair)
			except (IndexError, ValueError):
				break
	
	# Drop any duplicate entries (as a final safety measure)
	DF = DF.drop_duplicates(subset=['miRNA Name', 'Circ Name'], keep='first')
	logging.info(f"Final RH parsed entries: {len(DF)} (after deduplication)")
	DF.to_csv(hyb_parsed_file, index=False)
	
	if not keep_tmp:
		os.remove(hyb_tmp_file)

def TS_parser(fileIn, DB, acronym, keep_tmp=False):
    # Get the directory and base filename
    file_dir = os.path.dirname(fileIn)
    base_name = os.path.basename(fileIn)
    
    # Construct the input file path
    ts_tmp_file = os.path.join(file_dir, 'TS_TMP_' + base_name)
    
    # Skip if the file doesn't exist
    if not os.path.exists(ts_tmp_file):
        logging.warning(f"TargetScan output file {ts_tmp_file} not found. Skipping.")
        return DB
    
    # Read the file with pandas as in the old version
    df = pd.read_csv(ts_tmp_file, sep='\t')
    
    # Select and rename columns as in the old version
    df = df[['a_Gene_ID', 'miRNA_family_ID', 'UTR_start', 'UTR_end', 'Site_type']]
    df = df.rename(columns={
        'miRNA_family_ID': 'miRNA Name',
        'a_Gene_ID': 'Circ Name',
        'Site_type': 'Seed Category',
        'UTR_start': 'Start',
        'UTR_end': 'End'
    })
    
    # Prepend acronym to miRNA names as in the old version
    df['miRNA Name'] = acronym + '-' + df['miRNA Name']
    
    # Extract chromosome and strand from Circ Name
    df['Chrom'] = ''
    df['Strand'] = ''
    
    # Process each row to extract chromosome and strand
    for idx, row in df.iterrows():
        circ_name = row['Circ Name']
        try:
            # Format is typically chr:start-end or chr:start-end:strand
            circ_parts = circ_name.split(':')
            if len(circ_parts) >= 1:
                df.at[idx, 'Chrom'] = circ_parts[0]
            
            # Try to extract strand from the end
            if len(circ_parts) > 2 and circ_parts[-1] in ['+', '-']:
                df.at[idx, 'Strand'] = circ_parts[-1]
            else:
                # Try to extract from a pattern like chr1:12345-67890(+)
                strand_match = re.search(r'[()]([+-])[()]', circ_name)
                if strand_match:
                    df.at[idx, 'Strand'] = strand_match.group(1)
                else:
                    df.at[idx, 'Strand'] = '+'  # Default to + if not found
        except Exception as e:
            logging.warning(f"Error extracting chromosome/strand from {circ_name}: {e}")
    
    # Drop any duplicates before generating IDs
    df = df.drop_duplicates(subset=['miRNA Name', 'Circ Name'], keep='first')
    
    # Generate IDs as in the old version
    df['ID'] = np.arange(len(DB), (len(DB) + len(df)))
    df['ID'] = 'INT_TS_' + df['ID'].astype(str)
    
    # Cleanup temp file if needed
    if not keep_tmp:
        os.remove(ts_tmp_file)
    
    # Concatenate with DB
    DB = pd.concat([DB, df], ignore_index=True)
    
    return DB

def RH_parser(fileIn, DB=None, keep_tmp=False):
    # Get the directory and base filename
    file_dir = os.path.dirname(fileIn)
    base_name = os.path.basename(fileIn)
    
    # Construct the input file path
    hyb_parsed_file = os.path.join(file_dir, 'hyb_Parsed_' + base_name)

    # Create empty DataFrame with correct columns if DB is None
    if DB is None:
        DB = pd.DataFrame(columns=['Chrom', 'Start', 'End', 'miRNA Name', 'Circ Name', 'Strand', 'Seed Category', 'ID'])
    
    # Skip if the file doesn't exist
    if not os.path.exists(hyb_parsed_file):
        logging.warning(f"RNAhybrid parsed file {hyb_parsed_file} not found. Skipping.")
        return DB
    
    try:
        df = pd.read_csv(hyb_parsed_file)
        if not df.empty:
            # Drop any duplicates that might still exist after previous processing
            df = df.drop_duplicates(subset=['miRNA Name', 'Circ Name'], keep='first')
            
            # Extract chromosome and strand from Circ Name if available
            # Format is typically chr:start-end:strand
            if 'Circ Name' in df.columns:
                # Try to extract chromosome and strand from Circ Name
                try:
                    # Add Chrom column if not present
                    if 'Chrom' not in df.columns:
                        df['Chrom'] = df['Circ Name'].str.split(':', expand=True)[0]
                    
                    # Add Strand column if not present
                    if 'Strand' not in df.columns:
                        # Try to extract strand from the end of Circ Name
                        strand_extraction = df['Circ Name'].str.extract(r'([+-])$')
                        if not strand_extraction.empty and not strand_extraction.iloc[:, 0].isna().all():
                            df['Strand'] = strand_extraction.iloc[:, 0]
                        else:
                            # Default to '+' if strand can't be extracted
                            df['Strand'] = '+'
                except Exception as e:
                    logging.warning(f"Error extracting chromosome/strand from Circ Name: {e}")
                    # Set default values if extraction fails
                    if 'Chrom' not in df.columns:
                        df['Chrom'] = ''
                    if 'Strand' not in df.columns:
                        df['Strand'] = '+'
            
            # If DB is provided, prepare the ID column with appropriate offset
            if len(DB) > 0:
                start_idx = len(DB)
                df['ID'] = ['INT_RH_' + str(start_idx + i) for i in range(len(df))]
            else:
                # When used in parallel mode, we'll assign IDs later
                df['ID'] = ['INT_RH_temp_' + str(i) for i in range(len(df))]
            
            # Ensure all required columns are present
            required_columns = ['Chrom', 'Start', 'End', 'miRNA Name', 'Circ Name', 'Strand', 'Seed Category', 'ID']
            for col in required_columns:
                if col not in df.columns:
                    if col in ['Chrom', 'Strand']:
                        df[col] = ''  # Empty string for missing Chrom/Strand
                    elif col in ['Start', 'End']:
                        df[col] = 0   # Zero for missing Start/End
                    elif col == 'Seed Category':
                        df[col] = 'No Seed'  # Default seed category
                    else:
                        df[col] = 'Unknown'  # Default for other columns
            
            # Select columns in the right order, including all required columns
            new_df = df[required_columns]
            
            # Concatenate with DB
            result = pd.concat([DB, new_df], ignore_index=True)
        else:
            result = DB
        
        if not keep_tmp:
            os.remove(hyb_parsed_file)
    except (FileNotFoundError, pd.errors.EmptyDataError):
        result = DB

    return result

def single_coord(strand, coord_final, DB, circ):
	chunk = coord_final.loc[coord_final['name'] == circ]
	circ_chunk = DB.loc[DB['Circ Name'] == circ]
	if strand == '-':
		circ_chunk['fac_Start'] = chunk.iloc[0]['end'] - circ_chunk['End'].astype(int)
		circ_chunk['fac_End'] = chunk.iloc[0]['end'] - circ_chunk['Start'].astype(int)
	else:
		circ_chunk['fac_Start'] = circ_chunk['Start'].astype(int) + chunk.iloc[0]['start']
		circ_chunk['fac_End'] = circ_chunk['End'].astype(int) + chunk.iloc[0]['start']
	circ_chunk['Start'] = circ_chunk['fac_Start']
	circ_chunk['End'] = circ_chunk['fac_End']
	circ_chunk = circ_chunk.drop(['fac_Start','fac_End'], axis = 1)
	circ_chunk['Chrom'] = chunk.iloc[0]['chrom']
	circ_chunk['Strand'] = chunk.iloc[0]['strand']
	return circ_chunk

def multiple_coord(strand, coord_final, DB, circ):
	chunk = coord_final.loc[coord_final['name'] == circ]
	# Convert score to numeric but don't use it for coordinate calculations
	chunk['score'] = pd.to_numeric(chunk['score'])
	
	if strand == '-':
		chunk = chunk.reindex(index=chunk.index[::-1])
		circ_chunk = DB.loc[DB['Circ Name'] == circ].reset_index(drop=True)
		
		# First chunk: use end coordinate instead of score for filtering
		first_end = chunk['end'].iloc[0]
		# For negative strand, we need to find entries that fall within the first exon
		IDX = circ_chunk.index[circ_chunk['Start'].astype(float) < (first_end - chunk['start'].iloc[0])].to_list()
		df_c = circ_chunk.iloc[IDX, :]
		circ_chunk = circ_chunk.drop(circ_chunk.index[IDX]).reset_index(drop=True)
		
		# Calculate actual genomic coordinates for the first chunk
		df_c['fac_Start'] = chunk['end'].iloc[0] - df_c['End'].astype(int)
		df_c['fac_End'] = chunk['end'].iloc[0] - df_c['Start'].astype(int)
		df_c['Start'] = df_c['fac_Start']
		df_c['End'] = df_c['fac_End']
		df_c = df_c.drop(['fac_Start','fac_End'], axis = 1)
		df_c['Chrom'] = chunk.iloc[0]['chrom']
		df_c['Strand'] = chunk.iloc[0]['strand']
		
		# Process remaining chunks
		exon_length_so_far = chunk['end'].iloc[0] - chunk['start'].iloc[0]
		
		for i in range(1, len(chunk)):
			current_exon_length = chunk['end'].iloc[i] - chunk['start'].iloc[i]
			next_exon_boundary = exon_length_so_far + current_exon_length
			
			# Find entries that fall within this exon
			IDX = circ_chunk.index[circ_chunk['Start'].astype(float) < next_exon_boundary].to_list()
			df = circ_chunk.iloc[IDX, :]
			if df.empty:
				exon_length_so_far += current_exon_length
				continue
				
			circ_chunk = circ_chunk.drop(circ_chunk.index[IDX]).reset_index(drop=True)
			
			# Calculate actual genomic coordinates
			df['fac_Start'] = chunk['end'].iloc[i] - (df['End'].astype(int) - exon_length_so_far)
			df['fac_End'] = chunk['end'].iloc[i] - (df['Start'].astype(int) - exon_length_so_far)
			df['Start'] = df['fac_Start']
			df['End'] = df['fac_End']
			df = df.drop(['fac_Start','fac_End'], axis = 1)
			df['Chrom'] = chunk.iloc[i]['chrom']
			df['Strand'] = chunk.iloc[i]['strand']
			
			# Update running total of exon length
			exon_length_so_far += current_exon_length
			
			# Drop rows and columns with all NA values before concatenation
			df = df.dropna(how='all').dropna(axis=1, how='all')
			if not df.empty:
				df_c = pd.concat([df_c, df], ignore_index=True)
	else:
		circ_chunk = DB.loc[DB['Circ Name'] == circ].reset_index(drop=True)
		
		# First chunk: use start coordinate instead of score for filtering
		first_start = chunk['start'].iloc[0]
		first_end = chunk['end'].iloc[0]
		first_length = first_end - first_start
		
		# For positive strand, we need to find entries that fall within the first exon
		IDX = circ_chunk.index[circ_chunk['Start'].astype(float) < first_length].to_list()
		df_c = circ_chunk.iloc[IDX, :]
		circ_chunk = circ_chunk.drop(circ_chunk.index[IDX]).reset_index(drop=True)
		
		# Calculate actual genomic coordinates for the first chunk
		df_c['fac_Start'] = chunk['start'].iloc[0] + df_c['Start'].astype(int)
		df_c['fac_End'] = chunk['start'].iloc[0] + df_c['End'].astype(int)
		df_c['Start'] = df_c['fac_Start']
		df_c['End'] = df_c['fac_End']
		df_c = df_c.drop(['fac_Start','fac_End'], axis = 1)
		df_c['Chrom'] = chunk.iloc[0]['chrom']
		df_c['Strand'] = chunk.iloc[0]['strand']
		
		# Initialize df_c with proper columns if it's empty
		if df_c.empty:
			df_c = pd.DataFrame(columns=circ_chunk.columns)
		
		# Process remaining chunks
		exon_length_so_far = first_length
		
		# Process each chunk after the first one
		for i in range(1, len(chunk)):
			current_exon_length = chunk['end'].iloc[i] - chunk['start'].iloc[i]
			next_exon_boundary = exon_length_so_far + current_exon_length
			
			# Find entries that fall within this exon
			IDX = circ_chunk.index[circ_chunk['Start'].astype(float) < next_exon_boundary].to_list()
			df = circ_chunk.iloc[IDX, :]
			if df.empty:
				exon_length_so_far += current_exon_length
				continue
				
			circ_chunk = circ_chunk.drop(circ_chunk.index[IDX]).reset_index(drop=True)
			
			# Calculate actual genomic coordinates
			df['fac_Start'] = chunk['start'].iloc[i] + (df['Start'].astype(int) - exon_length_so_far)
			df['fac_End'] = chunk['start'].iloc[i] + (df['End'].astype(int) - exon_length_so_far)
			df['Start'] = df['fac_Start']
			df['End'] = df['fac_End']
			df = df.drop(['fac_Start','fac_End'], axis = 1)
			df['Chrom'] = chunk.iloc[i]['chrom']
			df['Strand'] = chunk.iloc[i]['strand']
			
			# Update running total of exon length
			exon_length_so_far += current_exon_length
			
			# Drop rows and columns with all NA values before concatenation
			df = df.dropna(how='all').dropna(axis=1, how='all')
			if not df.empty:
				df_c = pd.concat([df_c, df], ignore_index=True)
	
	return df_c

def add_interactions(db, interactions, column):
	if interactions.empty == False:
		for item in list(interactions[7].unique()):
			db.loc[db.ID == item, column] = 'Yes'
	return db

def calculate_exons(bedfile, gtf, rrna):
	logging.info("Loading input bed file: %s", str(bedfile))
	input_bed = pybedtools.BedTool(bedfile)
	bed_to_df = input_bed.to_dataframe()
	logging.info("Loading input GTF file: %s",str(gtf))
	gtf_file = pybedtools.BedTool(gtf)
	logging.info("Loading rRNA file: %s", str(rrna))
	rrna_file = pybedtools.BedTool(rrna)
	data = gtf_file.intersect(input_bed, wo=True, s=True)
	logging.info("Generating exons dataframe and polishing data")
	df_data = data.to_dataframe(names=['chr','gene_biotype','feature','start','end','6','strand','7','gene_info','c_chr','c_start','c_end','c_name','c_length','c_strand','overlap']).drop(['6','7','gene_biotype'], axis=1)
	exons = df_data.loc[df_data['feature'] == 'exon']
	exons = pd.concat([exons, exons['gene_info'].str.split(';', expand=True).rename(columns={0:'gene_id',1:'transcript_id',2:'exon_number',3:'gene_name',4:'gene_biotype',5:'transcript_name'})], axis=1).drop(['gene_info'], axis=1)
	exons['gene_id'] = exons['gene_id'].map(lambda x: x.lstrip(' gene_id ')).str.replace(r'"', '')
	exons['transcript_id'] = exons['transcript_id'].map(lambda x: x.lstrip(' transcript_id ')).str.replace(r'"', '')
	exons['exon_number'] = exons['exon_number'].map(lambda x: x.lstrip(' exon_number ')).str.replace(r'"', '')
	exons['gene_name'] = exons['gene_name'].map(lambda x: x.lstrip(' gene_name ')).str.replace(r'"', '')
	exons['gene_biotype'] = exons['gene_biotype'].map(lambda x: x.lstrip(' gene_biotype ')).str.replace(r'"', '')
	exons['transcript_name'] = exons['transcript_name'].map(lambda x: x.lstrip(' transcript_name ')).str.replace(r'"', '')
	exons['start'] = exons['start']-1
	logging.info("Returning dataframe")
	return exons, bed_to_df

def create_circular_DB(exons, bed_dataframe, rrna_file, genome, tmp_dir="./"):
    coordinates = pd.DataFrame()
    circulars = list(exons['c_name'].unique())
    logging.info("Calculating circulars coordinates")
    for c in circulars:
        trans = list(set(exons.loc[exons['c_name'] == c, 'transcript_id'].to_list()))
        for t in trans:
            slice = exons.loc[(exons['c_name'] == c)&(exons['transcript_id'] == t)].sort_values(by=['start'])
            if(slice['start'].iloc[0] == slice['c_start'].iloc[0]) and (slice['end'].iloc[-1] == slice['c_end'].iloc[-1]):
                coordinates = coordinates.append(slice[['chr','start','end','c_name','overlap','strand']])
    coordinates = coordinates.append(
        bed_dataframe.loc[
            bed_dataframe['name'].isin(
                list(set(list(bed_dataframe['name'].unique())).symmetric_difference(set(list(coordinates['c_name'].unique()))))
            )
        ].rename(columns={'chrom':'chr','name':'c_name','score':'overlap'})
    ).drop_duplicates()
    coord_bed = pybedtools.BedTool.from_dataframe(coordinates)
    coordinates = coord_bed.intersect(rrna_file, wa=True, v=True)
    coord_final = coordinates.to_dataframe()
    
    # Use a proper temporary file path
    fasta_output = os.path.join(tmp_dir, f'output_fasta_{os.getpid()}.fa')
    coordinates.sequence(fi=genome, fo=fasta_output, tab=True, name=True, s=True)
    fasted = pd.read_csv(fasta_output, sep='\t', names=['circs','sequence'])
    fasted['circs'] = fasted['circs'].str.replace(r"\(.*\)","",regex=True)
    try:
        os.remove(fasta_output)
    except OSError as e:
        logging.error(f"Failed to remove temporary file {fasta_output}: {e}")
    fasted = pd.concat([fasted, fasted['circs'].str.split('::', expand=True).rename(columns={0:'circ_name',1:'position'})], axis=1).drop(['circs'],axis=1)
    circulars = list(fasted['circ_name'].unique())
    logging.info("Retrieving FASTA sequences")
    Data = {}
    for circ_name in circulars:
        Data['>'+circ_name] = ''.join(fasted.loc[fasted['circ_name']==circ_name]['sequence'].to_list())
    return Data, circulars, coord_final

def starting_from_coord(bedfile, rrna, genome, tmp_dir="./"):
    logging.info("Loading input bed file with coordinates: %s", str(bedfile))
    input_bed = pybedtools.BedTool(bedfile)
    logging.info("Loading rRNA file: %s", str(rrna))
    rrna_file = pybedtools.BedTool(rrna)
    coordinates = input_bed.intersect(rrna_file, wa=True, v=True)
    coord_final = coordinates.to_dataframe()
    
    # Use a proper temporary file path
    fasta_output = os.path.join(tmp_dir, f'output_fasta_{os.getpid()}.fa')
    coordinates.sequence(fi=genome, fo=fasta_output, tab=True, name=True, s=True)
    fasted = pd.read_csv(fasta_output, sep='\t', names=['circs','sequence'])
    fasted['circs'] = fasted['circs'].str.replace(r"\(.*\)","",regex=True)
    try:
        os.remove(fasta_output)
    except OSError as e:
        logging.error(f"Failed to remove temporary file {fasta_output}: {e}")
    fasted = pd.concat([fasted, fasted['circs'].str.split('::', expand=True).rename(columns={0:'circ_name',1:'position'})], axis=1).drop(['circs'],axis=1)
    circulars = list(fasted['circ_name'].unique())
    logging.info("Retrieving FASTA sequences")
    Data = {}
    for circ_name in circulars:
        Data['>'+circ_name] = ''.join(fasted.loc[fasted['circ_name']==circ_name]['sequence'].to_list())
    return Data, circulars, coord_final

def getUTR(taxa):
	if taxa == 9606:
		return '3utr_human'
	elif taxa == 10090:
		return '3utr_human'
	elif taxa == 6239:
		return '3utr_worm'
	else:
		return '3utr_fly'

def process_file_with_parsers(output_file, tmp_dir, acronym, keep_tmp):
    """Process a single output file with all parsers in sequence."""
    try:
        file_path = os.path.join(tmp_dir, output_file)
        # Initialize an empty DataFrame for this file
        file_DB = pd.DataFrame(columns=['miRNA Name','Circ Name', 'Seed Category', 'Start', 'End', 'ID'])
        
        # Apply each parser in sequence
        file_DB = examine(file_path, file_DB, acronym, keep_tmp)
        file_DB = TS_parser(file_path, file_DB, acronym, keep_tmp)
        file_DB = RH_parser(file_path, file_DB, keep_tmp)
        return file_DB
    except FileNotFoundError:
        logging.warning(f'File {file_path} not found. Skipping')
        return pd.DataFrame(columns=['miRNA Name','Circ Name', 'Seed Category', 'Start', 'End', 'ID'])
    except Exception as e:
        logging.error(f"Error processing {output_file}: {str(e)}")
        return pd.DataFrame(columns=['miRNA Name','Circ Name', 'Seed Category', 'Start', 'End', 'ID'])


def third_party_processes(data_dict, taxa_number, cores, mature_file, acronym, miranda_sc, miranda_en, miranda_scale, miranda_go, miranda_ge, rnahybid_max, tmp_dir="./", keep_tmp=False, miR_family_info=None):
	#Set the processes for the run of the three software
	logging.info("Setting the processes for the three software run up")
	
	# Initialize the database DataFrame
	DB = pd.DataFrame(columns = ['miRNA Name','Circ Name', 'Seed Category', 'Start', 'End', 'ID'])
	
	processes = []
	out_list = []
	file_path = getCircrPath()
	UTR = getUTR(taxa_number)
	
	# Use provided miR_family_info or default
	if miR_family_info is None:
		miR_family_info = file_path + '/support_files/miRNA/miR_family_info.txt'
	
	# Ensure tmp_dir exists
	if not os.path.exists(tmp_dir):
		os.makedirs(tmp_dir)
	
	for item in data_dict.keys():
		output = str(item)[1:]+'.txt'
		out_list.append(output)
		inp = os.path.join(tmp_dir, str(item)[1:]+'_inputfile.fa')
		ts_inp = os.path.join(tmp_dir, 'TS_'+str(item)[1:]+'_inputfile.fa')
		
		# Define output files
		mir_output = os.path.join(tmp_dir, 'mir_TMP_'+output)
		ts_output = os.path.join(tmp_dir, 'TS_TMP_'+output)
		hyb_output = os.path.join(tmp_dir, 'hyb_TMP_'+output)
		hyb_parsed_file = os.path.join(tmp_dir, 'hyb_Parsed_'+output)
		
		# Only add processes if files don't already exist or keep_tmp is False
		if not keep_tmp or not (os.path.exists(mir_output) and os.path.exists(ts_output) and os.path.exists(hyb_output)):
			# Create input files only if we need to run the processes
			w = open(inp,'w')
			k = open(ts_inp,'w')
			w.write(item+'\n')
			w.write(data_dict[item]+'\n')
			k.write(item.lstrip('>')+'\t'+str(taxa_number)+'\t'+data_dict[item]+'\n')
			w.close()
			k.close()
			
			# Add processes to the list
			if not (keep_tmp and os.path.exists(mir_output)):
				processes.append('miranda {2} {0} -sc {3} -en {4} -scale {5} -go {6} -ge {7} -out {1}'.format(inp, mir_output, mature_file, miranda_sc, miranda_en, miranda_scale, miranda_go, miranda_ge))
			
			if not (keep_tmp and os.path.exists(ts_output)):
				processes.append('{0}/targetscan_70.pl {1} {2} {3}'.format(file_path, miR_family_info, ts_inp, ts_output))
			
			if not (keep_tmp and os.path.exists(hyb_output)):
				processes.append('RNAhybrid -m {4} -t {0} -q {2} -s {3} > {1}'.format(inp, hyb_output, mature_file, UTR, rnahybid_max))

	# Only run processes if there are processes to run
	if processes:
		#pool the analysis of the 3 software (RNAhybrid is the bottleneck) - only if there are more than one circular
		if len(out_list)==1:
			logging.info("Running miRNA:circRNA binding site predictions")
			for p in processes:
				run_process(p)
			
			# Only parse RNAhybrid output if the parsed file doesn't exist or keep_tmp is False
			if not (keep_tmp and os.path.exists(os.path.join(tmp_dir, 'hyb_Parsed_' + out_list[0]))):
				logging.info("Parsing RNAhybrid's output files")
				hybrid(os.path.join(tmp_dir, out_list[0]), keep_tmp)
			
			# Clean up input files if not keeping temporary files
			if not keep_tmp:
				for f in os.listdir(tmp_dir):
					if f.endswith('_inputfile.fa'):
						os.remove(os.path.join(tmp_dir, f))
		else:
			logging.info("Running miRNA:circRNA binding site predictions")
			pool = Pool(processes=cores)
			pool.map(run_process, processes)
			pool.close()
			pool.join()
			
			#clean input and generate new pool for RNAhybrid
			if not keep_tmp:
				for f in os.listdir(tmp_dir):
					if f.endswith('_inputfile.fa'):
						os.remove(os.path.join(tmp_dir, f))
			
			# Only parse hyb files that need parsing
			hyb_files_to_process = []
			for output in out_list:
				hyb_parsed_file = os.path.join(tmp_dir, 'hyb_Parsed_' + output)
				if not (keep_tmp and os.path.exists(hyb_parsed_file)):
					hyb_files_to_process.append(os.path.join(tmp_dir, output))
			
			if hyb_files_to_process:
				logging.info("Parsing RNAhybrid's output files")
				dpool = Pool(processes=cores)
				dpool.map(partial(hybrid, keep_tmp=keep_tmp), hyb_files_to_process)
				dpool.close()
				dpool.join()
			else:
				logging.info("Skipping RNAhybrid parsing as parsed files already exist")
	else:
		logging.info("Skipping miRNA:circRNA binding site predictions as all temporary files already exist")

	#create a single dataframe with all predicted circulars from 3 software
	logging.info("Collecting results into a single location")
	
	# Use multithreading to process files in parallel
	with ThreadPoolExecutor(max_workers=cores) as executor:
		# Submit all file processing tasks
		future_to_file = {
			executor.submit(process_file_with_parsers, output, tmp_dir, acronym, keep_tmp): output 
			for output in out_list
		}

		# Collect results as they complete
		all_dfs = []
		for future in concurrent.futures.as_completed(future_to_file):
			output = future_to_file[future]
			try:
				file_DB = future.result()
				if not file_DB.empty:
					all_dfs.append(file_DB)
			except Exception as e:
				logging.error(f"Error processing {output}: {str(e)}")

	# Combine all results and fix IDs
	if all_dfs:
		DB = pd.concat(all_dfs, ignore_index=True)
		
		# Fix the IDs to ensure they're sequential but preserve the original prefixes
		if not DB.empty and 'ID' in DB.columns:
			# Group by the software type in the ID
			id_parts = [id_val.split('_') for id_val in DB['ID']]
			
			# Create a dictionary to track counters for each prefix type
			prefix_counters = {}
			
			# Create new IDs that preserve the original prefix but use sequential numbers
			new_ids = []
			for i, parts in enumerate(id_parts):
				if len(parts) >= 2:
					# Extract the prefix (INT_TS, INT_M, INT_RH, etc.)
					prefix = f"{parts[0]}_{parts[1]}"
					
					# Initialize counter for this prefix if not seen before
					if prefix not in prefix_counters:
						prefix_counters[prefix] = 0
					
					# Create new ID with preserved prefix and sequential number
					new_id = f"{prefix}_{prefix_counters[prefix]}"
					prefix_counters[prefix] += 1
				else:
					# Fallback for malformed IDs
					new_id = f"INT_UNKNOWN_{i}"
				
				new_ids.append(new_id)
			
			# Update the IDs
			DB['ID'] = new_ids
	else:
		# Ensure DB is initialized with the correct columns
		DB = pd.DataFrame(columns=['miRNA Name','Circ Name', 'Seed Category', 'Start', 'End', 'ID'])
	
	#return the generated database
	return DB

def compare_predicted(coordinates, circulars, database):
	# Create an empty DataFrame with predefined columns
	tab = pd.DataFrame(columns = ['Chrom', 'Start', 'End', 'miRNA Name','Circ Name', 'Strand', 'Seed Category', 'ID'])
	for circ in circulars:
		result_df = None
		if (coordinates.loc[coordinates['name'] == circ, 'strand'].max() == '-') and (len(coordinates.loc[coordinates['name'] == circ]) == 1):
			result_df = single_coord('-', coordinates, database, circ)
		elif (coordinates.loc[coordinates['name'] == circ, 'strand'].max() == '-') and (len(coordinates.loc[coordinates['name'] == circ]) > 1):
			result_df = multiple_coord('-', coordinates, database, circ)
		elif (coordinates.loc[coordinates['name'] == circ, 'strand'].max() == '+') and (len(coordinates.loc[coordinates['name'] == circ]) == 1):
			result_df = single_coord('+', coordinates, database, circ)
		elif (coordinates.loc[coordinates['name'] == circ, 'strand'].max() == '+') and (len(coordinates.loc[coordinates['name'] == circ]) > 1):
			result_df = multiple_coord('+', coordinates, database, circ)
		
		if result_df is not None and not result_df.empty:
			# Ensure result_df has the same columns as tab
			result_df = result_df.reindex(columns=tab.columns)
			# Drop rows and columns with all NA values
			result_df = result_df.dropna(how='all').dropna(axis=1, how='all')
			if not result_df.empty:
				tab = pd.concat([tab, result_df], ignore_index=True)
	return tab

#Select miRanda's output when multiple softwares are able o find same interaction
#and report number of software that find that interaction
def process_duplicate_pair(duplicate_pair, tab):
    """Process a single duplicate pair in parallel with optimized performance."""
    miRNA_name, circ_name = duplicate_pair
    
    # Filter once to get all matching rows
    kkk = tab.loc[(tab['miRNA Name'] == miRNA_name) & (tab['Circ Name'] == circ_name)]
    
    if kkk.empty:
        return None
    
    # Extract software information once
    software_parts = [id_val.split('_') for id_val in kkk['ID']]
    software_types = [parts[1] if len(parts) > 1 else '' for parts in software_parts]
    
    # Get unique software types and count once
    unique_software = set(software_types)
    software_count = len(unique_software)
    
    # Make sure we preserve all columns, especially Chrom, Start, End, and Strand
    # We'll only exclude the columns we're going to split from the ID
    
    # Determine which row to use based on priority
    if 'M' in unique_software:
        # Find the first row with miRanda
        for idx, software in enumerate(software_types):
            if software == 'M':
                # Get the complete row to preserve all information
                result = kkk.iloc[idx].to_dict()
                result['Software Matched'] = software_count
                return result
    elif software_count > 1:
        # Find the row with highest priority software
        priority_software = min(unique_software)
        for idx, software in enumerate(software_types):
            if software == priority_software:
                # Get the complete row to preserve all information
                result = kkk.iloc[idx].to_dict()
                result['Software Matched'] = software_count
                return result
    elif software_count == 1 and len(kkk) > 1:
        # Just take the first row
        # Get the complete row to preserve all information
        result = kkk.iloc[0].to_dict()
        result['Software Matched'] = 1
        return result
    
    return None

def filter_interactions(valInt, AGO, tab, num_threads=None):
    """
    Filter interactions with parallel processing for better performance.
    
    Parameters:
    -----------
    valInt : str
        Path to validated interactions file
    AGO : str
        Path to AGO interactions file
    tab : pandas.DataFrame
        Table of interactions to filter
    num_threads : int, optional
        Number of threads to use for parallel processing. If None, uses all available cores.
    """
    import time
    start_time = time.time()
    
    # Determine number of threads to use
    if num_threads is None:
        num_threads = os.cpu_count()
    
    validated_interactions = pybedtools.BedTool(valInt)
    AGO_interactions = pybedtools.BedTool(AGO)
    
    # Add debug logging
    logging.info(f"Initial table shape: {tab.shape}")
    logging.info(f"First few rows of ID column: {tab['ID'].head().tolist()}")
    logging.info(f"Using {num_threads} threads for parallel processing")
    
    # Process duplicates as before
    filter_start = time.time()
    filtered = tab.merge(tab[tab.duplicated(subset=['miRNA Name', 'Circ Name'], keep=False)][['miRNA Name', 'Circ Name']], how='left', indicator=True)
    filtered = filtered[filtered['_merge'] == 'left_only'].drop(['_merge'],axis=1)
    filtered['Software Matched'] = 1
    filter_end = time.time()
    logging.info(f"Initial filtering took {filter_end - filter_start:.2f} seconds")
    
    # Add debug logging
    logging.info(f"Number of non-duplicate entries: {filtered.shape[0]}")
    
    # Process duplicates and build to_add list as before
    duplicates_start = time.time()
    duplicates = list(set(list(tab[tab.duplicated(subset=['miRNA Name', 'Circ Name'], keep=False)][['miRNA Name', 'Circ Name']].apply(tuple, axis=1))))
    duplicates_end = time.time()
    logging.info(f"Finding duplicates took {duplicates_end - duplicates_start:.2f} seconds")
    
    # Add debug logging
    logging.info(f"Number of duplicate pairs found: {len(duplicates)}")
    
    # Process duplicates in parallel
    to_add = []
    if duplicates:
        parallel_start = time.time()
        with ThreadPoolExecutor(max_workers=num_threads) as executor:
            # Submit all duplicate processing tasks
            future_to_duplicate = {
                executor.submit(process_duplicate_pair, duplicate, tab): duplicate 
                for duplicate in duplicates
            }
            
            # Collect results as they complete
            for future in concurrent.futures.as_completed(future_to_duplicate):
                duplicate = future_to_duplicate[future]
                try:
                    result = future.result()
                    if result is not None:
                        to_add.append(result)
                except Exception as e:
                    logging.error(f"Error processing duplicate {duplicate}: {str(e)}")
        parallel_end = time.time()
        logging.info(f"Parallel processing of {len(duplicates)} duplicates took {parallel_end - parallel_start:.2f} seconds")
        logging.info(f"Average time per duplicate: {(parallel_end - parallel_start) / len(duplicates):.4f} seconds")
    
    # Add the processed items to filtered
    concat_start = time.time()
    if len(to_add) > 0:
        logging.info(f"Adding {len(to_add)} processed duplicate entries")
        to_add_df = pd.DataFrame(to_add)
        # Drop rows and columns with all NA values before concatenation
        to_add_df = to_add_df.dropna(how='all').dropna(axis=1, how='all')
        if not to_add_df.empty:
            filtered = pd.concat([filtered, to_add_df], ignore_index=True)
    else:
        logging.info("No processed duplicates to add")
    concat_end = time.time()
    logging.info(f"Concatenating results took {concat_end - concat_start:.2f} seconds")
    
    # Ensure proper types for BED format
    # Convert Start and End columns to integers
    filtered['Start'] = filtered['Start'].astype(int)
    filtered['End'] = filtered['End'].astype(int)
    
    # Ensure we have the essential BED columns in the right order
    bed_columns = ['Chrom', 'Start', 'End', 'miRNA Name', 'Circ Name', 'Strand']
    for col in bed_columns:
        if col not in filtered.columns:
            logging.error(f"Required column {col} not found in DataFrame")
            # Create a minimal valid BED file to avoid errors
            return filtered
    
    try:
        # Create BedTool from the DataFrame with proper column order
        bed_start = time.time()
        filtered_bed = pybedtools.BedTool.from_dataframe(filtered[bed_columns + [c for c in filtered.columns if c not in bed_columns]])
        bed_end = time.time()
        logging.info(f"Creating BedTool took {bed_end - bed_start:.2f} seconds")
        
        # Define a function to process intersections in parallel
        def process_intersection(intersection_type, bed_file, filtered_bed, filtered):
            intersection_start = time.time()
            if not os.path.exists(bed_file) or os.path.getsize(bed_file) == 0:
                logging.warning(f"{intersection_type} file {bed_file} is empty or doesn't exist. Skipping intersection.")
                return None, None
            
            try:
                # Create BedTool for the intersection file
                intersection_bed = pybedtools.BedTool(bed_file)
                
                # Perform the intersection
                intersected = filtered_bed.intersect(intersection_bed, wa=True, wb=True, f=0.50)
                intersected_df = intersected.to_dataframe(names=list(range(0, len(filtered.columns) + len(intersection_bed.to_dataframe().columns))))
                
                # Process intersection results
                intersected_df = intersected_df.iloc[:, 0:len(filtered.columns)]
                intersected_df.columns = filtered.columns
                
                intersection_end = time.time()
                logging.info(f"{intersection_type} intersection took {intersection_end - intersection_start:.2f} seconds")
                return intersection_type, intersected_df
            except Exception as e:
                logging.error(f"Error during {intersection_type} intersection: {e}")
                return None, None
        
        # Process intersections in parallel
        intersection_results = {}
        intersect_start = time.time()
        with ThreadPoolExecutor(max_workers=2) as executor:  # Only 2 intersections, so max 2 threads
            # Submit intersection tasks
            future_to_intersection = {
                executor.submit(process_intersection, "Validated", valInt, filtered_bed, filtered): "Validated",
                executor.submit(process_intersection, "AGO", AGO, filtered_bed, filtered): "AGO"
            }
            
            # Collect results as they complete
            for future in concurrent.futures.as_completed(future_to_intersection):
                intersection_type = future_to_intersection[future]
                try:
                    result_type, result_df = future.result()
                    if result_type is not None and result_df is not None:
                        intersection_results[result_type] = result_df
                except Exception as e:
                    logging.error(f"Error processing {intersection_type} intersection: {str(e)}")
        intersect_end = time.time()
        logging.info(f"Parallel intersections took {intersect_end - intersect_start:.2f} seconds")
        
        # Process validated interactions result
        merge_start = time.time()
        if "Validated" in intersection_results:
            val_intersected = intersection_results["Validated"]
            val_intersected['Validated'] = 'Yes'
            filtered = filtered.merge(val_intersected[['miRNA Name', 'Circ Name', 'Validated']], how='left', on=['miRNA Name', 'Circ Name'])
            filtered['Validated'] = filtered['Validated'].fillna('No')
        else:
            filtered['Validated'] = 'No'
        
        # Process AGO interactions result
        if "AGO" in intersection_results:
            AGO_intersected = intersection_results["AGO"]
            AGO_intersected['AGO'] = 'Yes'
            filtered = filtered.merge(AGO_intersected[['miRNA Name', 'Circ Name', 'AGO']], how='left', on=['miRNA Name', 'Circ Name'])
            filtered['AGO'] = filtered['AGO'].fillna('No')
        else:
            filtered['AGO'] = 'No'
        merge_end = time.time()
        logging.info(f"Merging intersection results took {merge_end - merge_start:.2f} seconds")
    
    except Exception as e:
        logging.error(f"Error creating BedTool from DataFrame: {e}")
        # Add the Validated and AGO columns with default 'No' values
        filtered['Validated'] = 'No'
        filtered['AGO'] = 'No'
        
    # Final check for any remaining duplicates
    dedup_start = time.time()
    logging.info(f"Before final deduplication: {filtered.shape[0]} rows")
    # Use all columns except ID for deduplication check
    id_cols = [col for col in filtered.columns if col != 'ID']
    filtered = filtered.drop_duplicates(subset=id_cols, keep='first')
    dedup_end = time.time()
    logging.info(f"After final deduplication: {filtered.shape[0]} rows")
    logging.info(f"Final deduplication took {dedup_end - dedup_start:.2f} seconds")
    
    end_time = time.time()
    total_time = end_time - start_time
    logging.info(f"Total filter_interactions processing time: {total_time:.2f} seconds")
    
    return filtered

#Add the notation from the CircBase database dump file
#from the supprot files
def add_CircBase_annotation(table, infile, genome_version, circbase_annot=None):
	#get path of the CircBase ref file
	file_path = getCircrPath()
	
	# Use provided circbase_annot or default
	if circbase_annot is None:
		circbase_annot = file_path + '/support_files/circBase_circRNA.txt'
	
	CircBase = pd.read_csv(circbase_annot, sep='\t', names=['Chrom', 'Start', 'End', 'CircName', 'Strand', 'Organism'])
	input_bed = pd.read_csv(infile, sep='\t', names=['Chrom', 'Start', 'End', 'CircName', 'Score', 'Strand'])
	#selecting only circulars of the same genome version investigated
	CircBase = CircBase.loc[CircBase['Organism'] == genome_version]
	#creating the reference column
	CircBase['ToBeMatched'] = CircBase['Chrom'] + CircBase['Start'].astype(str) + CircBase['End'].astype(str) + CircBase['Strand'].astype(str)
	# creating the reference dictionary from circBase
	CircBase_dict = dict(zip(CircBase.ToBeMatched, CircBase.CircName))
	match_dict = {}
	for ix, xg in input_bed.groupby(['Chrom', 'CircName', 'Strand']):
		xg = xg.sort_values(['Start'])
		start = xg['Start'].iloc[0]
		end = xg['End'].iloc[-1]
		# creating the string to be matched in circBase
		tbm = ix[0] + start.astype(str) + end.astype(str) + ix[2]
		# if there is a match in circBase
		if tbm in CircBase_dict:
			# add circBase ID
			match_dict[ix[1]] = CircBase_dict[tbm]
		else:
			# else, add 'Not Available'
			match_dict[ix[1]] = 'NA'

	# now update the output interaction file
	table['circBase ID'] = table['Circ Name']
	table['circBase ID'] = table['circBase ID'].map(match_dict)
	return table

def fix_bed_scores(input_bed, output_bed=None):
    """
    Fix the score column in a BED file by calculating the length of each interval
    when the score is missing or invalid.
    
    Parameters:
    -----------
    input_bed : str
        Path to the input BED file
    output_bed : str, optional
        Path to save the processed BED file. If None, will modify the input file in-place.
        
    Returns:
    --------
    str
        Path to the processed BED file
    """
    logger = logging.getLogger(__name__)
    logger.info(f"Fixing scores in BED file: {input_bed}")
    
    # Read the BED file into a pandas DataFrame
    try:
        # Try to read with standard BED columns
        bed_df = pd.read_csv(input_bed, sep='\t', header=None)
    except Exception as e:
        logger.error(f"Error reading BED file: {e}")
        return input_bed
    
    # Ensure we have at least 3 columns (chrom, start, end)
    if bed_df.shape[1] < 3:
        logger.error("BED file must have at least 3 columns")
        return input_bed
    
    # If we have at least 5 columns (score is column 4, 0-indexed)
    if bed_df.shape[1] >= 5:
        # Convert score column to numeric, with errors='coerce' to convert non-numeric to NaN
        bed_df[4] = pd.to_numeric(bed_df[4], errors='coerce')
        
        # For rows with NaN scores, calculate the length of the interval
        mask = bed_df[4].isna()
        if mask.any():
            logger.info(f"Found {mask.sum()} rows with invalid scores")
            
            # Calculate length (end - start) for rows with invalid scores
            bed_df.loc[mask, 4] = bed_df.loc[mask, 2] - bed_df.loc[mask, 1]
            
            logger.info("Replaced invalid scores with interval lengths")
    
    # If no output file is specified, create a temporary file
    if output_bed is None:
        temp_file = tempfile.NamedTemporaryFile(delete=False, suffix='.bed')
        output_bed = temp_file.name
        temp_file.close()
    
    # Write the processed DataFrame to the output file
    bed_df.to_csv(output_bed, sep='\t', header=False, index=False)
    
    # If we created a temporary file, replace the original with it
    if output_bed != input_bed and output_bed != temp_file.name:
        logger.info(f"Saving processed BED file to: {output_bed}")
    elif output_bed == temp_file.name:
        shutil.move(output_bed, input_bed)
        output_bed = input_bed
        logger.info(f"Updated original BED file: {input_bed}")
    
    return output_bed

def main(bedfile, gtf, rrna, AGO, TAXA, outfile, Genome, cores, valid_interactions, mirna_mature, mirna_acronym, coord, genome_version,
		 miranda_sc, miranda_en, miranda_scale, miranda_go, miranda_ge, rnahybrid_max, tmp_dir="./", keep_tmp=False, miR_family_info=None, circbase_annot=None):  # valid_interactions
	begin = time.time()
	temp_files = []  # Track temporary files for cleanup
	
	try:
		logging.info("Starting Circr analysis pipeline")
		
		# Fix scores in the BED file
		logging.info("Fixing BED scores in input file")
		bedfile = fix_bed_scores(bedfile)
		
		if coord:
			logging.info("Creating Circular Database and Retrieving Coordinates from provided coordinates")
			circ_dict, circulars, coordinates = starting_from_coord(bedfile, rrna, Genome, tmp_dir)
		else:
			logging.info("Retrieving exon information from GTF file")
			data_exons, df_bed = calculate_exons(bedfile, gtf, rrna)
			logging.info("Creating Circular Database and Retrieving Coordinates from exon data")
			circ_dict, circulars, coordinates = create_circular_DB(data_exons, df_bed, rrna, Genome, tmp_dir)
		
		logging.info("Performing analysis with RNAhybrid, miRanda and TargetScan")
		circular_database = third_party_processes(circ_dict, TAXA, cores, mirna_mature, mirna_acronym, 
			miranda_sc, miranda_en, miranda_scale, miranda_go, miranda_ge, rnahybrid_max, tmp_dir, keep_tmp, miR_family_info)
		#print("circular_database", circular_database)
		logging.info("Merging overlapping seed predictions")
		table = compare_predicted(coordinates, circulars, circular_database)
		#print("table", table)
		logging.info("Adding Annotation to Predicted Interactions")
		interactions = filter_interactions(valid_interactions, AGO, table, num_threads=cores)
		#print("interactions", interactions)
		logging.info('Adding information from CircBase annotation')
		interactions = add_CircBase_annotation(interactions, bedfile, genome_version, circbase_annot)
		
		logging.info('Saving results table to %s', outfile)
		interactions.to_csv(outfile, index=False)
		
		end = time.time()
		logging.info('Total execution time: %s seconds', str(end - begin))
		logging.info('Analysis complete. Closing Circr')
		logging.info('Circr analysis performed in %s minutes', str(round((end - begin) / 60)))
		
		return interactions
		
	except FileNotFoundError as e:
		logging.error(f"File not found error: {e}")
		logging.error(f"Please check that all input files exist and are accessible")
		traceback.print_exc()
		end = time.time()
		logging.info('Execution time before error: %s seconds', str(end - begin))
		raise
	except pd.errors.EmptyDataError as e:
		logging.error(f"Empty data error: {e}")
		logging.error(f"One of the input files appears to be empty")
		traceback.print_exc()
		end = time.time()
		logging.info('Execution time before error: %s seconds', str(end - begin))
		raise
	except MemoryError:
		logging.error("Memory error: Not enough memory to complete the analysis")
		logging.error("Try using a machine with more RAM or reducing the size of your input files")
		end = time.time()
		logging.info('Execution time before error: %s seconds', str(end - begin))
		raise
	except Exception as e:
		logging.error(f"Unexpected error during Circr analysis: {str(e)}")
		traceback.print_exc()
		end = time.time()
		logging.info('Execution time before error: %s seconds', str(end - begin))
		raise
	finally:
		# Clean up any temporary files if not keeping them
		if not keep_tmp:
			for temp_file in temp_files:
				if os.path.exists(temp_file):
					try:
						os.remove(temp_file)
						logging.debug(f"Removed temporary file: {temp_file}")
					except Exception as e:
						logging.warning(f"Failed to remove temporary file {temp_file}: {e}")

if __name__ == "__main__":
	######## PARSER SECTION ###########
	parser = argparse.ArgumentParser(
	    prog='python3 Circr.py',
	    formatter_class=argparse.RawDescriptionHelpFormatter,
	    description="Circr help section")
	parser.add_argument("-i","--input", help="List of Circular RNA or their exon coordinates.")
	parser.add_argument("-c","--coord", action="store_true", help="Defines if the input file contains exon coordinates rather than Circular RNA coordinates")
	parser.add_argument("-s","--organism", type=str,help="Defines the selected organism for analysis. Available organisms are: human, mouse, worm, fruitfly. Default is human.", default="human")
	parser.add_argument("-v","--genome_version", type=str, help="Defines the genome version of the selected organism. Versions available are: human (hg19, hg38), mouse (mm9, mm10), worm (ce10, ce11), fruitfly (dm3, dm6). Default for human is hg38", default="none")
	parser.add_argument("--gtf", type=str,help="Alternative location for GTF file for the analysis. Default uses data for selected organism.", default="none")
	parser.add_argument("--genome", type=str,help="Alternative location for genome file for the analysis. Default uses data for selected organism.", default="none")
	parser.add_argument("--rRNA", type=str,help="Alternative location for ribosomal RNA file for the analysis. Default uses data for selected organism.", default="none")
	parser.add_argument("--miRNA", type=str,help="Alternative location for miRNA file for the analysis. Default uses data for selected organism.", default="none")
	parser.add_argument("--miR_family_info", type=str,help="Alternative location for miRNA family information file. Default uses data for selected organism.", default="none")
	parser.add_argument("--circbase_annot", type=str,help="Alternative location for CircBase annotation file. Default uses data for selected organism.", default="none")
	parser.add_argument("--AGO", type=str, help="Alternative AGO peaks file. Default uses data for selected organism.", default="none")
	parser.add_argument("--validated_interactions", type=str, help="Alternative validated interactions file. Default uses data for selected organism.", default="none")
	parser.add_argument("--threads", type=int, help="Set the number of threads for multiprocess. Default is 8 cores", default=8)
	parser.add_argument("-o","--output", type=str,help="Defines output file name. Default is Circr_Analysis.csv", default="Circr_Analysis.csv")
	parser.add_argument("-Msc", type=float, help="Set score threshold, default 140. miRanda Parameter", default="140.0")
	parser.add_argument("-Men", type=float, help="Set energy threshold to -value Kcal/mol. Default 1. miRanda Parameter", default="1.0")
	parser.add_argument("-Mscale", type=float, help="Set scaling parameter. Default 4. miRanda Parameter", default="4.0")
	parser.add_argument("-Mgo", type=float, help="Set gap-open penalty to -value. Default 4 (-4). miRanda Parameter", default="-4.0")
	parser.add_argument("-Mge", type=float, help="Set gap-extend penalty to -value. Default 9 (-9). miRanda Parameter", default="-9.0")
	parser.add_argument("-RHmax", type=int, help="The  maximum  allowed  length of a target sequence. Default is 10000000. RNAhybrid Parameter", default="10000000")
	parser.add_argument("--tmp_dir", type=str, help="Directory for temporary files. Default is current directory.", default="./")
	parser.add_argument("--log", type=str, help="Log file path. Default is run.log", default="run.log")
	parser.add_argument("--keep_tmp", action="store_true", help="Keep temporary files after analysis.")

	args = parser.parse_args()

	logging.basicConfig(filename=args.log, format='%(asctime)s - %(levelname)s:	%(message)s', datefmt='%m/%d/%Y %I:%M:%S %p', level=logging.INFO)

	### VARIABLE CHECKS ###
	if (args.genome_version == 'none') and (args.organism == 'human'):
		args.genome_version = 'hg38'

	if args.organism not in organisms.keys():
		logging.error("Organism not supported. Exiting...")
		sys.exit()

	if args.genome_version not in organisms[args.organism]:
		logging.error("Genome version not supported. Exiting...")
		sys.exit()

	### SETTING DEFAULTS ###

	if args.gtf == 'none':
		args.gtf = getCircrPath() + DataTree[args.organism][args.genome_version]['GTF']
	if args.genome == 'none':
		args.genome = getCircrPath() + DataTree[args.organism][args.genome_version]['Genome']
	if args.rRNA == 'none':
		args.rRNA = getCircrPath() + DataTree[args.organism][args.genome_version]['rRNA_Coords']
	if args.miRNA == 'none':
		args.miRNA = getCircrPath() + DataTree[args.organism]['file']
	if args.miR_family_info == 'none':
		args.miR_family_info = getCircrPath() + '/support_files/miRNA/miR_family_info.txt'
	if args.circbase_annot == 'none':
		args.circbase_annot = getCircrPath() + '/support_files/circBase_circRNA.txt'
	if args.AGO == 'none':
		args.AGO = getCircrPath() + DataTree[args.organism][args.genome_version]['AGO']
	if args.validated_interactions == 'none':
		args.validated_interactions = getCircrPath() + DataTree[args.organism][args.genome_version]['INT']

	try:
		targetscan_path = f"{getCircrPath()}/targetscan_70.pl"
		subprocess.run(['chmod', '+x', targetscan_path], check=True)
		logging.info(f"Successfully set permissions for targetscan_70.pl")
	except subprocess.CalledProcessError as e:
		logging.error(f"Failed to set permissions for targetscan_70.pl: {e}")
		sys.exit(1)

	logging.info('Running analysis on %s; genome set to %s and version %s', args.input, args.organism, args.genome_version)

	main(args.input,
		 args.gtf,
		 args.rRNA,
		 args.AGO,
		 DataTree[args.organism][args.genome_version]['TAXA'],
		 args.output,
		 args.genome,
		 args.threads,
		 args.validated_interactions,
		 args.miRNA, DataTree[args.organism]['acronym'],
		 args.coord,
		 args.genome_version,
		 args.Msc,
		 args.Men,
		 args.Mscale,
		 args.Mgo,
		 args.Mge,
		 args.RHmax,
		 args.tmp_dir,
		 args.keep_tmp,
		 args.miR_family_info,
		 args.circbase_annot)