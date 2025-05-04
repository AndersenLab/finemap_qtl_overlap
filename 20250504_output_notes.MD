Aim: Add output statements to the main.nf workflow using updated NF output methods

Add statement before the workflow definition to publish outputs
```
// Needed to publish results
nextflow.preview.output = true
```

Within the workflow block, we need to add a publish statement with this syntax
```
publish:
{process_id}.{outputname}
```
To publish the output of the `calculate_ld` process `ld_files`

```
publish:
  calculate_ld.ld_files
```

Then after the workflow defintion this block
```
output {
    "." {
        mode "copy"
    }
}
```
