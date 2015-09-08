#!/bin/bash   


time for c in {1..22}
do
	for p in {0..5}
	do
		echo "chromosome "$c", population "$p
		python mutation_sets.py -pop $p -c $c
	done
done
