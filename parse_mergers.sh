# for line in $(cat ./checks/mergers_1.txt)

while read line
do
    if [[ "$(echo $line | cut -c1-1)" = "#" ]]; then
        echo "skipping first line"
        continue
    fi
    final_scale=$(echo $line | awk '{print $1}')
    final_sat_mass=$(echo $line | awk '{print $2}')
    final_cen_mass=$(echo $line | awk '{print $3}')
    final_sat_id=$(echo $line | awk '{print $4}')
    final_cen_id=$(echo $line | awk '{print $5}')

    sat_id=$final_sat_id
    cen_id=$final_cen_id
    echo "\n\nnew satellite"
    echo $sat_id
    echo $cen_id
    while [ $sat_id -ne $cen_id ]
    do
        grep_output_sat=$(grep "^ *[0-1]\.[0-9]* *$sat_id " tree_0_0_0.dat)
        grep_output_cen=$(grep "^ *[0-1]\.[0-9]* *$cen_id " tree_0_0_0.dat)

        echo $(echo $grep_output_sat | awk '{ print $1, $2, $10, $11, $15, $3, $4}')
        echo $(echo $grep_output_cen | awk '{ print $1, $2, $10, $11, $15, $3, $4}')

        sat_id=$(echo $grep_output_sat | awk '{print $4}')
        cen_id=$(echo $grep_output_cen | awk '{print $4}')

        echo 
    done
    # 

done < ./checks/mergers_1.txt