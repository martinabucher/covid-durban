match_score=1.000000
mismatch_score= -1.000000

target_internal_open_gap_score= -0.8500000
target_internal_extend_gap_score= -0.500000
target_left_open_gap_score= 0.000000
target_left_extend_gap_score= 0.000000
target_right_open_gap_score= 0.000000
target_right_extend_gap_score= 0.000000
query_internal_open_gap_score= -1.3
query_internal_extend_gap_score= -0.9
query_left_open_gap_score= -1.1
query_left_extend_gap_score= -0.8
query_right_open_gap_score= -0.45
query_right_extend_gap_score=-0.25

query_internal_open_gap_score= -1.3
query_internal_extend_gap_score= -0.9
query_left_open_gap_score= -1.1
query_left_extend_gap_score= -0.8
query_right_open_gap_score=   query_left_open_gap_score 
query_right_extend_gap_score= query_left_extend_gap_score

target_internal_open_gap_score  =query_internal_open_gap_score
target_internal_extend_gap_score=query_internal_extend_gap_score
target_left_open_gap_score=      0.0
target_left_extend_gap_score=    0.0
target_right_open_gap_score=     target_left_open_gap_score
target_right_extend_gap_score=   target_left_extend_gap_score

def alignerInit(aligner):
  aligner.match_score=match_score
  aligner.mismatch_score=mismatch_score
  aligner.target_internal_open_gap_score=target_internal_open_gap_score
  aligner.target_internal_extend_gap_score=target_internal_extend_gap_score
  aligner.target_left_open_gap_score=target_left_open_gap_score
  aligner.target_left_extend_gap_score =target_left_extend_gap_score
  aligner.target_right_open_gap_score= target_right_open_gap_score 
  aligner.target_right_extend_gap_score= target_right_extend_gap_score
  aligner.query_internal_open_gap_score= query_internal_open_gap_score
  aligner.query_internal_extend_gap_score= query_internal_extend_gap_score
  aligner.query_left_open_gap_score= query_left_open_gap_score
  aligner.query_left_extend_gap_score= query_left_extend_gap_score
  aligner.query_right_open_gap_score=query_right_open_gap_score
  aligner.query_right_extend_gap_score=query_right_extend_gap_score
  return 

