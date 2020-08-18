#!/bin/bash
# Usage: remove all utility bills pdf file password 
shopt -s nullglob
declare -a arr=("AL" "AK" "AZ" "AR" "CA" "CO" "CT" "DE" "DC" "FL" "GA"  "HI" "ID" "IL" "IN" "IA" "KS" "KY" "LA" "ME" "MD" "MA" "MI" "MN" "MS" "MO" "MT" "NE" "NV" "NH" "NJ" "NM" "NY" "NC" "ND" "OH" "OK" "OR" "PA" "RI" "SC" "SD" "TN" "TX" "UT" "VT" "VA" "WA" "WV" "WI" "WY")
for i in "${arr[@]}"
do
   echo "working on state $i"
   Rscript USA_UR_Deaths_LancetID.R $i
done

