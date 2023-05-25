## Splicekit and the singularity container

There are several software dependencies for splicekit, and we packaged them into a single singularity container.

If you would like to build the container from the provided [splicekit.def](splicekit.def) file you need to have `singularity` installed. Please follow the operating system specific [official singularity documentation](https://docs.sylabs.io/guides/3.0/user-guide/installation.html).

## Building the splicekit singularity image

The [splicekit.def](splicekit.def) file contains all the third-party open-source software splicekit calls to do the analysis.

To build the splicekit singularity image run:

`$ singularity build splicekit.sif splicekit.def`

This will create the `splicekit.sif` container file. When you obtain the `splicekit.sif` file, you can relocate it to whichever folder you like.

Finally, in the `splicekit.config` file, you should modify the `container` variable to `singularity run path_to/splicekit.sif`.
