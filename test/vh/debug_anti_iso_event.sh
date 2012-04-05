#!/bin/bash

# Run megadebug on all of the UW events
# If you save this to a file, view it with less -r

megadebug=$fsa/TMegaSelector/scripts/megadebug 

echo "Debugging MMMT"

$megadebug AnalyzeMMMT.py zh_events_ntuple.txt /mmmt/final/Ntuple \
  unique- base_selections+  hadronic_tau_id not_m3_id \
  --events 176201,271,421507855 

# all
  #171050,643,821417992 \
  #179434,371,641824424 \
  #179889,78,81028181 \
  #179889,168,235516917 \
  #175975,360,515584039 \
  #176201,271,421507855 \
  #176467,146,226304277 \
  #177718,788,1263752726 \
  #177730,189,264203514 \
  #178098,32,42470279 \
  #166512,1063,1231425135 \
  #166841,62,66008312 \
  #167281,272,346878062 \
  #167786,107,126622509 
