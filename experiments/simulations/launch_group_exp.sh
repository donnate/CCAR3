#!/bin/bash

# Define the values for the variables
#!/bin/bash

# Define the values for the variables
theta_strengths="high medium low"
#in_values="500"
r_values="3"
#p_values="100 300 500 700"
p_values="100 300 500 700 1000 2000"
q_values="10 30 50 80"
for theta in $theta_strengths; do
  for p in $p_values; do
      for r in $r_values; do
	  for q in $q_values; do
             sbatch rrr_group_experiment.sh "$theta" "$r" "$p" "$q"
          done
    done
  done
done
