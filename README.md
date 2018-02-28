## mendacity

simply-specified tests for inheritance models.

## syntax

All tests are specified in `YAML` format.

Each test file specifies one pedigree and a number of tests for that pedigree.

Depths, GQ's, PL's, etc are not specified only the number of alternates:

| genotype | n-alts | name    |
| -------- | ------ | ------- |
|  0/0     |   0    | HOM-REF |
|  0/.     |   0    | HOM-REF |
|  ./0     |   0    | HOM-REF |
|  0/1     |   1    | HET     |
|  ./1     |   1    | HET     |
|  1/.     |   1    | HET     |
|  1/1     |   2    | HOM-ALT |
|  ./.     |  -1    | UNKNOWN |

For decomposed alleles `./.` can also be 0-alts.

A simple test for a trio of `mom`, `dad`, `proband` looks like:

```yaml
  name: ar-test
  description: simple autosomal-recessive
  alts: [[1, 1, 2]] 
  modes: ["autosomal-recessive"]
  not-modes: ["autosomal-dominant"]
```

where:
+ `modes` defines what inheritance modes this meets
+ `not-modes` defines what inheritance modes this must not meet

For this trio, the yaml specifying the pedigree is:

```yaml
pedigree:
- &dad
    id: Fred
    status: unaffected
    sex: male
- &mom
    id: Alice
    status: unaffected
    sex: female
-
    id: Bobby
    status: affected
    sex: female
    father: *dad
    mother: *mom
```

and the test-cases are:

```yaml
cases:
-
  name: ar-test
  description: simple autosomal-recessive
  alts: [[1, 1, 2]]
  modes: ["autosomal-recessive"]
  not-modes: ["autosomal-dominant"]
-
  name: comp-het-example
  description: >
     requires 2 alts from different parents
     and alleles.                                                                                                                                                                      
  modes: ["compound-het"]
  alts: [[0, 1, 1], [1, 0, 1]]
```

## constraints

+ status ∈ unaffacted|affected|unknown|carrier
+ sex ∈ male|female|unknown
+ Modes ∈ autosomal-dominant|autosomal-recessive|compound-het|de-novo|x-linked-recessive|x-linked-dominant|x-linked-de-novo

## VCF/PED

A VCF and PED pair can be created from a test-case yaml with:

```
python mendacity/tovcf.py tests/transmitted-dn.yaml
```

This will create `transmitted-dn.vcf` and `transmitted-dn.ped`



