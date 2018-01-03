#!/bin/bash

for OUTPUT in $(ls Chips_images/*eps)
do
	echo $OUTPUT
  convert $OUTPUT $OUTPUT.png
  rm $OUTPUT
done
