#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 21 09:46:36 2020

@author: thomas
"""

import sys
import copy
import urllib.request
import re
from statistics import mode, stdev

from collections import defaultdict
from Bio import SeqFeature, SeqRecord, AlignIO
from Bio.Align import AlignInfo

from Bio import BiopythonWarning
import warnings
with warnings.catch_warnings():
	warnings.simplefilter('ignore', BiopythonWarning)

def parse_alignment(path, frame):
	#path, frame = [args.match_alignment, args.force_alignment_frame]
	frame = frame - 1 if frame is not None else None
	
	alignment = AlignIO.read(path, "fasta")
	
	# Get modal start and finish positions
	prelim = dict()
	for seq_record in alignment:
		prelim[seq_record.id] = {e: len(re.search(r, str(seq_record.seq)).group()) for e, r in zip(["start", "finish"],["^-*", "-*$"])}
	start, finish = zip(*[[d[e] for e in ["start", "finish"]] for d in prelim.values()])
		
	modes = dict()
	for e, v in zip(['start', 'finish'], [start, finish]):
		modes[e] = mode(v) if frame is None else round(mode(v)/3)*3+frame
	
	# Generate consensus
	alignment_summary = AlignInfo.SummaryInfo(alignment)
	consensus = str(alignment_summary.gap_consensus())
	
	# Find approximate ungapped distance to modal position for each sequence
	output = dict()
	
	for seq_record in alignment:
		dists = dict()
		for e in ['start', 'finish']:
			#seq_record = alignment[298]
			#e = 'start'
			#e = 'finish'
			# Set up the two locations in a list
			locations = [prelim[seq_record.id][e], modes[e]]
			
			if(locations[0] == locations[1]):
				dists[e] = 0
			else:
				# Extract the sequence, either the target if it's outside the alignment or the consensus if inside
				sequence = str(seq_record.seq) if locations[0] < locations[1] else consensus
				
				# Reverse the sequence if we're looking at the finish
				sequence = sequence[::-1] if e == 'finish' else sequence
				
				# Extract the in between sequence
				between_sequence = sequence[min(locations):max(locations)]
				
				# Find the ungapped distance
				distance = len(between_sequence) - between_sequence.count("-")
				
				# Correct the distance depending on direction and end
				distance = distance * -1 if(locations[0] > locations[1]) else distance
				distance = distance * -1 if e == 'finish' else distance
				
				dists[e] = distance
				
			
		
		output[seq_record.id] = dists
	
	# Find standard deviation of the distance
	start, finish = zip(*[[d[e] for e in ["start", "finish"]] for d in output.values()])
	
	std_dev = {e:round(stdev(v)) for e, v in zip(["start", "finish"], [start, finish])}
	
	return(output, std_dev)

def str_is_int(s):
	try: 
		int(s)
		return True
	except ValueError:
		return False

def parse_overlap(overlap_in, annotation_in):
	sys.stderr.write("Running overlap autocorrection for %s based on:\n"  % (annotation_in))
	output = {}
	
	for overlapstring in overlap_in:
		#overlapstring = args.overlap[0]
		overlap_list = overlapstring.split(",")
		err = "Error: overlap string " + overlapstring
		if(len(overlap_list) != 2):
			sys.exit(err + " does not contain exactly two comma-separated elements")
		elif(not str_is_int(overlap_list[1])):
			sys.exit(err + " does not have an integer as the second item")
		sys.stderr.write("\toverlap of %s bp with %s\n" % tuple(reversed(overlap_list)))
		output[overlap_list[0]] = int(overlap_list[1])
	
	return(output)

def parse_stringsearch(annotation_in, startstring_in, finishstring_in, translation_table, match_alignment):
	#annotation_in, startstring_in, finishstring_in, translation_table, match_alignment = [args.annotation[0], args.startstring, args.finishstring, args.translation_table, args.match_alignment is not None]
	output = dict()
	sys.stderr.write("Running string searching for %s using:\n" % (annotation_in))
	
	for end, searchstring in zip(['start','finish'], [startstring_in, finishstring_in]):
		if(searchstring is None):
			continue
		if(end is 'start' and translation_table is None):
			sys.exit("Error: --translation_table must be specified if searching for --startstring")
		sys.stderr.write("\t%s: %s\n" % (end, searchstring))
		
		search = tuple(searchstring.split(","))
		
		err = "Error: search string " + searchstring
		if(search[0] not in ['A', 'N']):
			sys.exit(err + " has an unrecognised sequence type")
		elif(search[0] == 'A'):
			if((match_alignment and len(search) !=2) or (not match_alignment and len(search) != 3)):
				sys.exit(err + " does not have the correct number of items")
			elif(any(s not in list("/GPAVLIMCFYWHKRQNEDST*") for s in search[1])):
				sys.exit(err + " is specified as amino acid but non-standard character(s) included (must be GPAVLIMCFYWHKRQNEDST*)")
			elif(not match_alignment and str_is_int(search[2])):
				sys.exit(err + " is searching for an amino acid sequence but includes what seems to be a reading frame")
			elif(translation_table is None):
				sys.exit(err + " is specified as amino acid but no --translation_table given")
		else:
			if((match_alignment and len(search) !=3) or (not match_alignment and len(search) != 4)):
				sys.exit(err + " does not have the correct number of items")
			elif(any(s not in list("/ATGC") for s in search[1])):
				sys.exit(err + " is specified as nucleotide but non-standard character(s) included (must be ATGC)")
			elif(search[2] not in ['1','2','3','*']):
				sys.exit(err + " has an unrecognised reading frame")
		
		if(not match_alignment):
			mmhi = search[3] if search[0] == 'N' else search[2]
			if( mmhi not in ['F','L','C','FC','LC','N'] ):
				sys.exit(err + " has an unrecognised multiple-hit handling instruction ( must be F, FC, C, LC, L or N)")
		else:
			search = search + tuple(['C'])
		
		output[end] = search
	
	return(output)

def find_all(a_str, sub):
	start = 0
	while True:
		start = a_str.find(sub, start)
		if start == -1: return
		yield start
		start += 1

def loadnamevariants():
	output = {}
	for line in urllib.request.urlopen("https://raw.githubusercontent.com/tjcreedy/biotools/master/gene_name_variants.txt"):
		line = line.decode('utf-8').strip()
		name = line.split(";")[0]
		annotype = line.split(":")[0].split(";")[1]
		variants = line.split(":")[1].split(",")
		for v in variants:
			for g in ['', ' ']:
				v = v.replace(g, '')
				for s in ['',' GENE', ' '+annotype.upper()]:
					output[v+s] = name
	return(output)

def get_features_from_names(seqrecord, names, namevariants):
	#seqrecord, names = [seq_record, feature_names]
	
	names = [names] if isinstance(names, str) else names
	
	features = defaultdict(list)
	unrecognised_names = set()
	unidentifiable_features = set()
	#seqname = seqrecord.name
	
	for feat in seqrecord.features:
		#feat = seqrecord.features[68]
		# Remove any translations
		if('translation' in feat.qualifiers.keys()):
			del(feat.qualifiers['translation'])
		
		# Extract the feature name
		featname = 0
		nametags = ['gene', 'product', 'label', 'standard_name']
		if(any(t in feat.qualifiers.keys() for t in nametags)):
			for t in nametags:
				if(t in feat.qualifiers.keys()):
					featname = feat.qualifiers[t][0].upper()
					break
		elif(feat.type in ['source', 'misc_feature', 'repeat_region', 'D-loop', 'rep_origin','gap']):
			continue
		else:
			unidentifiable_features.add((feat.type, feat.location.start, feat.location.end))
			#err = "Warning, can't identify %s %s annotation" % (seqname, feat.type)
			#if(hasattr(feat.location, 'start')):
			#	err += " %s-%s" % (str(int(feat.location.start)+1), str(int(feat.location.end)))
			#err += ": no gene/product/label/standard_name tag\n"
			#sys.stderr.write(err)
			continue
		
		if(featname in namevariants):
			name = namevariants[featname]
			if(name in names):
				features[name].append(feat)
		else:
			unrecognised_names.add(featname)
	
	return(features, unrecognised_names, unidentifiable_features)

def syncronise_features(name, feats, synctype, seqname):
	#synctype = 'gene'
	#name, feats = list(features.items())[2]
	# Organise features
	target_feats = []
	other_feats = []
	
	for feat in feats:
		if(feat.type == synctype):
			target_feats.append(feat)
		else:
			other_feats.append(feat)
	
	warnstart = "Warning, " + seqname + " does not have any "
	warnend = " annotations of " + name + "\n"
	if(len(target_feats) < 1):
		#sys.stderr.write(warnstart + synctype + warnend)
		return
	if(len(other_feats) < 1):
		#sys.stderr.write(warnstart + "non-" + synctype + warnend)
		return
	
	# Check that the other features are all the same
	if(not all(other_feats[0].location == feat.location for feat in other_feats)):
		sys.stderr.write("Warning, positions of " + str(len(other_feats)) + " non-" + synctype + " annotations for " + name + " in " + seqname + " do not match, this entry will not be modified\n")
		return
		
	# Correct the target features
	for feat in target_feats:
		feat.location = other_feats[0].location

def correct_positions_by_overlap(target, context_features, overlap, maxoverlap, seqlength, seqname):
	#target, maxoverlap, seqlength = [feat, args.overlap_maxdist, len(seq_record.seq)]
	
	context_overdist = set()
	outfeat = copy.deepcopy(target)
	corrected_start, corrected_finish = [target.location.start, target.location.end]
	
	# Work through combinations of target and context features
	
	tpos = [int(target.location.start), int(target.location.end)]
	for context_name in context_features:
		#context_name = list(context_features.keys())[0]
		context = context_features[context_name][0]
		cpos = [int(context.location.start), int(context.location.end)]
		
		# Find gap between context and target
		gap = [min(cpos)-max(tpos), min(tpos)-max(cpos)]
		gap = [g for g in gap if g > 0]
		
		# Check gap is permissible
		if(len(gap) > 0 and gap[0] > maxoverlap + overlap[context_name]):
			context_overdist.add(context_name)
			continue
		
		# Set orientation (+ve, context follows target)
		orientation = 0
		# Set current distance between selected positions
		distance = 0
		# Set the index of the target position to change
		target_tpos_i = None
		
		
		# Find cross-wise distances
		crossdist = [max(cpos) - min(tpos), min(cpos) - max(tpos)]
		
		# Find the structure of the overlap
		if(abs(crossdist[0]) == abs(crossdist[1])):
			sys.stderr.write("Warning: orientation of " + context_name + " and " + "target annotations completely match, cannot perform overlap correction\n")
			continue
		elif(abs(crossdist[1]) < abs(crossdist[0])):
			# Target is before context, overlap should be latter position of target and first position of context
			orientation = 1
			distance = crossdist[1]
			target_tpos_i = tpos.index(max(tpos))
		else:
			# Target is after context, overlap should be first position of target and latter position of context
			orientation = -1
			distance = crossdist[0]
			target_tpos_i = tpos.index(min(tpos))
		
		# Calculate the exact new position
		corrected_tpos = tpos[target_tpos_i] + distance + orientation * (overlap[context_name])
		
		# Ensure the position is not outside the contig
		corrected_tpos = 0 if corrected_tpos < 0 else corrected_tpos
		corrected_tpos = seqlength if corrected_tpos > seqlength else corrected_tpos
		
		# Check that the orientation hasn't been flipped
		if((target_tpos_i == 0 and corrected_tpos >= int(target.location.end)) or (target_tpos_i == 1 and corrected_tpos <= int(target.location.start))):
			continue
		
		# Overwrite the relevant target end position
		if(target_tpos_i == 0):
			corrected_start = corrected_tpos
		else:
			corrected_finish = corrected_tpos
	
	outfeat.location = SeqFeature.FeatureLocation(corrected_start, corrected_finish, outfeat.location.strand)
	
	return(outfeat, context_overdist)

def correct_feature_by_alignment(feat, query_spec, distances):
	#query_spec, distances = [stringspec, alignment_distances[seqname]]
	
	outfeat = copy.deepcopy(feat)
	
	for end in ['start', 'finish']:
		#end = 'start'
		#end = 'finish'
		
		# Skip if not processing this end or if the end is already correct
		if(end not in query_spec.keys() or distances[end] == 0):
			continue
		
		if((end == 'start' and feat.location.strand == 1) or 
		   (end == "finish" and feat.location.strand == -1)):
			
			if(str_is_int(str(feat.location.start))): # Check if exact position, if not skip
				outfeat.location = SeqFeature.FeatureLocation(outfeat.location.start + distances[end], outfeat.location.end, outfeat.location.strand)
			
		else:
			
			if(str_is_int(str(feat.location.end))):
				outfeat.location = SeqFeature.FeatureLocation(outfeat.location.start, outfeat.location.end + distances[end], outfeat.location.strand)
	
	return(outfeat)

def correct_feature_by_query(feat, query_spec, seq_record, seqname, distance, featurename, translation_table, prioritise_longer_finishes):
	#query_spec, distance, featurename, translation_table, prioritise_longer_finishes = [stringspec, args.search_distance, name, args.translation_table, True]
	#feat, query_spec, distance, featurename, translation_table, prioritise_longer_finishes = [currfeat, stringspec, args.search_distance, name, args.translation_table, False]
	feat_start, feat_finish = feat.location.start, feat.location.end
	errstart = "Warning: sequence " + seqname + " has "
	
	for end in ['start','finish']:
		#end = 'start'
		#end = 'finish'
		
		if(end not in query_spec.keys()):
			continue
		
		# Unpack search tuple
		code, query, out_rf, selector = query_spec[end] if len(query_spec[end]) == 4 else query_spec[end] + tuple("X")
		selector = out_rf if code == 'A' else selector
		
		# Check if already ends with the searched sequence - REMOVED AS EXISTING SEQUENCE MAY BE SHORTER SUBSET OF DESIRED SEQUENCE e.g. TA TAA
		#if(end_already_correct(feat.extract(seq_record.seq), query, end, code, out_rf, args.translation_table)):
		#	continue
		
		# Generate sequence for searching
		subject_sequence, subject_start = extract_subject_region(seq_record, feat, end, code, distance)
		
		# Find locations of query
		results = dict()
		for q in query.split("/"):
			# Work through locations of hits
			for i in list(find_all(str(subject_sequence), q)):
				# If the current hit location exists and is longer than the current hit, do not change, else add current hit length
				results[i] = results[i] if i in results.keys() and len(q) < results[i] else len(q)
		
		# Remove results if selector is N, FC or LC
		errmid = ""
		if(len(results) == 0):
			errmid = "no matches of "
		if(selector == "N" and len(results) > 1):
			errmid = "multiple matches of "
			results = {}
		if(selector == "FC"):
			results = {i:l for i, l in results.items() if i <= distance}
			errmid = "no first closest matches of "
		elif(selector == "LC"):
			results = {i:l for i, l in results.items() if i >= distance}
			errmid = "no last closest matches of "
		
		# Retain only locations in the specified reading frame
		if(code == 'N' and out_rf != "*" and len(results) > 0):
			if(end == "start"):
				results = {i:l for i, l in results.items() if (i + distance) % 3 + 1 == int(out_rf)}
				# Retain location if that location's rf (l+1)%3 is equal to the (target rf converted to subject rf)
			else:
				results = {i:l for i, l in results.items() if (abs(feat.location.end - feat.location.start) - distance + i - 1) % 3 + 1 == int(out_rf)}
			errmid = "no matches in the specified frame of " if len(results) == 0 else errmid
		
		truncated = False
		codon_start = 1 if end == "start" else None
		
		if(end == "start"):
			
			# First check if truncated
			start_distance = len(seq_record) - feat.location.end if feat.location.strand == -1 else feat.location.start
			truncated = start_distance < distance
			
			# Retain only start locations that generate realistic amino acid sequences
			stopcounts = [get_stopcount(i, l, feat.location.strand, code, end, distance, subject_start, feat_start, feat_finish, seq_record, translation_table) for i, l in results.items()]
			results = {i:l for (i, l), c in zip(results.items(), stopcounts) if c <= 1 and c == min(stopcounts)}
			
			# If no feasible results at the start position, instead find the closest in-frame position 
			if(len(results) == 0):
				errmid = "no ORF-producing matches (will set to closest ORF) of "
				
				# Set the current position - if normal, this is the current start position, if truncated, this is the end of the contig
				contig_start = list(find_all(str(subject_sequence), 'N'))[-1] + 1 if truncated else None
				current_position = contig_start if truncated else distance
				
				# Set the correction - if normal, this is 0, if truncated, this is 1, to ensure no results outside the contig
				correction = 1 if truncated else 0
				
				# Set up the three alternative results
				results = { current_position + v + correction : 1 for v in [-1, 0, 1] }
				
				# Find the result with suitable ORF
				stopcounts = [get_stopcount(i, l, feat.location.strand, code, end, distance, subject_start, feat_start, feat_finish, seq_record, translation_table) for i, l in results.items()]
				results = {i:l for (i, l), c in zip(results.items(), stopcounts) if c <= 1 and c == min(stopcounts)}
				
				# If truncated, set result to contig start but note codon position
				if(len(results) > 0 and truncated):
						codon_start = sorted(results.keys())[0] - contig_start + 1
						results = {contig_start : 1}
				
			
		else:
			# Prioritise longer matches by removing any matches shorter than the longest match
			if(len(results) > 0 and prioritise_longer_finishes):
				max_length = max(results.values())
				results = {i:l for i, l in results.items() if l == max_length}
			
			
			# If searching for finish string and the annotation is likely truncated, remove any incomplete stop codons
			
			# Find the distance to the finish of the contig from the current position (strand-dependent)
			finish_distance = feat.location.start if feat.location.strand == -1 else len(seq_record) - feat.location.end
			
			# Remove partial stops if annotation likely truncated
			if(finish_distance < distance):
				truncated = True
				results = {i:l for i,l in results.items() if l > 2}
				# If no results remain, set result to be the end of the contig
				if(len(results) == 0):
					errmid = "no >2 matches near truncated finish (will set to contig end) of "
					results = {str(subject_sequence).find('N') - 1 : 1}
		
		# Parse location results
		errend = query  + " at the " + end + " of " + featurename + "\n"
		
		if(len(results) > 0):
			
			# Select the first location by default (which also is the LC location if LC set)
			locations = sorted(results.keys())
			location = locations[0]
			
			if(len(locations) > 1):
				if(selector in ['L', 'FC']):
					# If more than 1 locations but the user wants the last or first closest selected
					location = locations[-1]
				elif(selector == 'C'):
					loc_dist = [abs(distance-l) for l in locations]
					if(loc_dist.count(min(loc_dist)) == 1):
						location = locations[loc_dist.index(min(loc_dist))]
					else:
						errmid = "multiple closest matches (taking first) of "
			
			feat_start, feat_finish = get_newends(location, results[location], feat.location.strand, end, distance, code, subject_start, feat_start, feat_finish, truncated)
		else:
			errmid = "no succesful matches of " if errmid == "" else errmid
		
		if(errmid != ""):
			sys.stderr.write(errstart + errmid + errend)
	
	outfeat = copy.deepcopy(feat)
	outfeat.location = SeqFeature.FeatureLocation(feat_start, feat_finish, feat.location.strand)
	
	return(outfeat, codon_start)

def get_newends(location, length, strand, end, distance, code, subject_start, feat_start, feat_finish, truncated):
	
	# Convert location if on reverse strand
	location = location if strand == 1 else abs(location - (2*distance + 1))
	
	# Generate the new end position
		# Correct by length of the match if at the finish end
	change = location + strand * length if end == "finish" else location
		# Multiply by 3 if AA
	change = change * 3 if code == 'A' else change
		# Calculate
	newend = subject_start + change
	
	# Apply new end to appropriate end
	if((end == "start" and strand == 1) or (
			end == "finish" and strand == -1)):
		feat_start = SeqFeature.BeforePosition(newend) if truncated else SeqFeature.ExactPosition(newend)
	else:
		feat_finish = SeqFeature.AfterPosition(newend) if truncated else SeqFeature.ExactPosition(newend)
	
	return(feat_start, feat_finish)

def stopcount(seq_record, table, frame = (1,2,3)):
	
	# Check input types
	run_frame = (frame,) if not isinstance(frame, (tuple, list)) else frame
	
	# Run counting
	with warnings.catch_warnings():
		warnings.simplefilter('ignore', BiopythonWarning)
		counts = [seq_record.seq[(i-1):].translate(table = table).count("*") for i in run_frame]
	
	# Return string or list depending on length
	if(len(counts) > 1):
		return counts
	else:
		return counts[0]

def get_stopcount(location, length, strand, code, end, distance, subject_start, feat_start, feat_finish, seq_record, table):
	#location, length, strand, table = [19, results[20], feat.location.strand, args.translation_table]
	# Generate the potential end position for this location
	potstart, potfinish = get_newends(location, length, strand, end, distance, code, subject_start, feat_start, feat_finish, False)
	
	# Build the potential new feature for this location
	potfeat = SeqFeature.SeqFeature(SeqFeature.FeatureLocation(potstart, potfinish), strand = strand)
	
	# Extract the sequence for this potential new feature
	potseq = SeqRecord.SeqRecord(potfeat.extract(seq_record.seq))
	
	# Count the stops in this sequence and return true if less than or equal to 1
	return(stopcount(potseq, table, 1))



def extract_subject_region(seqrecord, feat, end, code, distance):
	'''For a given feature and end ("start" or "finish"), extract a sequence of either nucleotides (code = 'N') or amino acids (code = 'A'), consisting of the first or last position plus or minus positions equal to distance in the reading direction of the feature'''
	#seqrecord = seq_record
	
	# Find the centre point and distances for the subject region
	featend = "featstart"
	central_position = feat.location.start
	distances = (distance, distance + 1)
	if((end == "start" and feat.location.strand == -1) or (end == "finish" and feat.location.strand == 1)):
		featend = "featfinish"
		central_position = feat.location.end
		distances = (distance + 1, distance)
	
	# Convert for amino acids
	if(code == 'A'):
		# Multiply by 3
		distances = tuple(3*d for d in distances)
		
		#Find how many trailing out-of-frame bases
		trailing_bases = len(feat.extract(seqrecord.seq)) % 3
		
		if(end == "finish"):
			modifier = 0
			if(trailing_bases != 0):
				modifier = 3 - trailing_bases
			# Correct distances to ensure in-frame translation for out-of-frame trailing bases and for strand direction
			distances = (distances[0] + (-1 * feat.location.strand * modifier), distances[1] + (feat.location.strand * modifier))
	
	# Delimit the region
	end_positions = (central_position - distances[0], central_position + distances[1])
	
	# Truncate the region if it exceeds the contig
	start_position = 0 if end_positions[0] < 0 else end_positions[0]
	finish_position = len(seqrecord) if end_positions[1] > len(seqrecord) else end_positions[1]
	
	# Truncate the region if it exceeds the other end of the annotation
	if(featend == "featstart"):
		finish_position = feat.location.end if finish_position > feat.location.end else finish_position
	else:
		start_position = feat.location.start if start_position < feat.location.start else start_position
	
	# Generate the feature	
	subject_feat = SeqFeature.SeqFeature(SeqFeature.FeatureLocation(start_position, finish_position), strand = feat.location.strand)

	# Extract the sequence
	sequence = subject_feat.extract(seqrecord.seq)
	correction = [start_position - end_positions[0], end_positions[1] - finish_position]
	correction = correction[::-1] if feat.location.strand == -1 else correction
	sequence = correction[0] * "N" + sequence + correction[1] * "N"
	
	sequence = sequence.translate(table = 5) if(code == 'A') else sequence
	
	return(sequence, end_positions[0])

def correct_truncated_features(feature, contiglength):
	#feature, contiglength = [feat, len(seq_record)]
	
	corrected_start, corrected_finish = [feature.location.start, feature.location.end]
	
	if(int(feature.location.start) == 0):
		corrected_start = SeqFeature.BeforePosition(corrected_start)
	
	if(int(feature.location.end) == contiglength):
		corrected_finish = SeqFeature.AfterPosition(corrected_finish)
	
	feature.location = SeqFeature.FeatureLocation(corrected_start, corrected_finish, feature.location.strand)
	

def end_already_correct(nuc_seq, query_seq, end, code, frame, translation_table):
	#nuc_seq = feat.extract(seq_record.seq)
	#query_seq = query
	#frame = out_rf
	
	
	if(code == 'N'):
		return(any(end == "start" and frame == "1" and nuc_seq.startswith(q) or # Starts with sequence, rf is 1
			   (end == "finish" and nuc_seq.endswith(q) and nuc_seq.rfind(q) % 3 + 1 == int(frame) )) for q in query_seq.split("/"))
	else:
		aa_seq = nuc_seq.translate(table = translation_table)
		return(any((end == "start" and aa_seq.startswith(q)) or
		           (end == "finish" and aa_seq.endswith(query_seq)) for q in query_seq.split("/")))

def check_context_features(context_features, seqname):
	# Check context_features all match in positions
	for name, feats in context_features.items():
		locations = [feat.location for feat in feats]
		
		if(not all(locations[0] == loc for loc in locations)):
			
			err = "Warning, positions of " + str(len(locations)) + " annotations for " + str(name) + " in " + str(seqname) + " do not match. If these are multiple distinct loci, this should be fine, but if these should cover the same locus, you may get incorrect results.\n"
			
			for i, feat in enumerate(feats):
				err += "\t(" + str(i+1) + ") " + feat.type +" is located at bases " + str(int(feat.location.start)+1) + " to " + str(int(feat.location.end)) + "\n"
			
			sys.stderr.write(err)