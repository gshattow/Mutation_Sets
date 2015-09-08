#!/bin/bash   


time for c in {1..22}
do
	echo "parsing individuals for chromosome "$c
	./individuals.bash -c $c
done
