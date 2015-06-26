rna seq workflow for use with qproject.

Add section "params" in the config that is created by qproject like this:

```
{
    "base": "project",
    "data": "project/data",
    "ref": "project/ref/",
    "src": "project/src",
    "var": "project/var",
    "result": "project/result",
    "run": "project/run",
    "etc": "project/etc",
    "logs": "project/logs",
    "archive": "project/archive",
    "usr": "project/usr",
    "params": {
        "gtf": "genes.gtf",
        "indexedGenome": "Bowtie2Index/genome"
    }
}
```

where `indexedGenome` and `gtf` are paths relative to `ref`.

`indexedGenome` is the basename of a bowtie2 index.
