# GFFPal

A quick and flexible toolset/python library for working with GFF files.

I'm often needing to convert formats into or from GFF3.

Genometools and Aegean are useful tools for some aspects of GFF manipulation (especially tidying and conversion to GTF).
But they often aren't great at handling non-specifcation compliant GFF3s (i.e. CDS attached directly to genes, "match" style GFFs, or augustus hints GFFs).
Also I wanted something that I could use to select by ID (instead of tricky `awk` or `grep -f` commands).

This is intended for my own use and I'm motivated by ease of use rather than computational effeciency (hence Python).
But if you have thoughts, contributions, or whatever, please raise an issue to get in touch.
