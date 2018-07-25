import os, re
import subprocess
import shlex
import sys


cmd="convert -quality 60 /Users/taweewat/Documents/red_sequence/pisco_color_plots/redsq_8.55_Field292_tilted.png /Users/taweewat/Documents/red_sequence/pisco_color_plots/star_galaxy_sep_12_2Field292.png /Users/taweewat/Documents/red_sequence/pisco_color_plots/psf_est/psf_est_Field292.png output.pdf"
print cmd
sub = subprocess.check_call(shlex.split(cmd))
