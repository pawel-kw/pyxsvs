#! /bin/bash
rm -rf *.png
rm -rf *.edf
rm -rf SiD_1-1_40C_report.pdf
rm -rf *.p
sed --in-place '/mask\.edf/d' xsvs_input.txt 
