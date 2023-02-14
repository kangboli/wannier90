../wannier90.x aiida
tail -n +2 tdc.debug.print_bk_perm_test | awk '{print $3-$7,$4-$8,$5-$9}'| uniq
