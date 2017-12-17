# UA HPC Workflow

## ICE

Automatically submit single core jobs until simulation goal is complete.

`checkque_ice.sh`

Run `checkque_ice.sh` as
```bash
/home/u15/agladstein/ABC/macsSwig_AJmodels/checkque_ice.sh GOAL QUEMAX VERSION MODEL SYSTEM
```
Where `GOAL` is the number of simulations you want to complete,  
`QUEMAX` is the number of jobs allowed in the que,  
`VERSION` is the type of simulation (e.g. `instant`, `genome`),  
`MODEL` is the model number (`1`, `2`, or `3`)
`SYSTEM` is the ICE cluster (`cluster`, `smp`, `htc`)

for example:
```bash
/home/u15/agladstein/ABC/macsSwig_AJmodels/HPC_workflow/checkque_ice.sh 500000 200 genome 1 cluster
```

In order to run the switch files must be present in the `HPC_workflow` directory,
```bash
/home/u15/agladstein/ABC/macsSwig_AJmodels/HPC_workflow/switch1.txt
/home/u15/agladstein/ABC/macsSwig_AJmodels/HPC_workflow/switch2.txt
/home/u15/agladstein/ABC/macsSwig_AJmodels/HPC_workflow/switch3.txt
```

`checkque_ice.sh` -> `main_function_AJmodel_j2.sh` -> `chain_pbs.j2`

### Template PBS

`main_function_AJmodel_j2.sh` uses jinja templates to create the appropriate PBS scripts

There are four jinja templates:  
`chain_pbs.j2`  
`macsargs_pbs.j2`  
`run-sim_pbs.j2`  
`genome_stats_pbs.j2`


## Ocelote