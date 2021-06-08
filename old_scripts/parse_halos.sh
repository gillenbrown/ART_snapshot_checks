for central_id in $(awk {'print $3'} ./checks/growth_1.txt)
do
    grep "^ *[0-1]\.[0-9]* *$central_id " tree_0_0_0.dat | awk '{ print $1, $2, $10, $11, $15, $3, $4}'
done
