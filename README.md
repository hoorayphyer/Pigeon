# Pigeon

When electrons and protons meet a rotating [neutron star](https://en.wikipedia.org/wiki/Neutron_star) possessing a magnetic field literally millions of millions of [Earth's own magnetic field](https://en.wikipedia.org/wiki/Earth%27s_magnetic_field), what does it look like?

<iframe src="https://drive.google.com/file/d/1Sx84V2R3XWCLOb7f9JLNRw1kei9rm0EL/preview" width="640" height="480" allow="autoplay"></iframe>
<br>(simulation of a pulsar, modeled as a rapidly rotating neutron star surrounded with charged particles. The simulation was done on 1120 cores)

The technique to scientifically (as well as artistically) portrait the above is [Particle-in-Cell](https://en.wikipedia.org/wiki/Particle-in-cell) simulation, or PIC simulation. It is widely used in studies of astrophysical and laboratory plasma phenomena. Pigeon is yet another implementation of this technique. Although born out of a specific research project on [pulsars](https://en.wikipedia.org/wiki/Pulsar), Pigeon is meant to be a general-purpose PIC simulator (whatever that is), and embraces the following:
* object-oriented programming
* open source
* reproducible research
* multi-tasking

## Technical specifics
* supported OS: Linux, OS X
* written in progressively modern-er C++
* data saved with [HDF5](https://en.wikipedia.org/wiki/Hierarchical_Data_Format)
* TDD with ample unit tests using [Catch2](https://github.com/catchorg/Catch2)
* fully parallelized using [MPI](https://en.wikipedia.org/wiki/Message_Passing_Interface), and has successfully run using 120
octacosa-core(28-core) Intel Broadwell processors on the Pleiades
supercomputer at NASA’s Ames Research Center
* achieves virtually perfect load balancing by overlaying the traditional
Cartesian topology with a [primary/replica](https://en.wikipedia.org/wiki/Master/slave_(technology)) paradigm
* supports simulation resumption
* offers high flexibility of customization for specific physics problems
* supports runtime configuration via a TOML file, parsed by [`toml++`](https://github.com/marzer/tomlplusplus/).

## Labs and shadows
Have ten projects at hand? Don't download a Pigeon for each, create a lab for each instead! 
Powered by CMake, Pigeon supports "one codebase multiple builds". A lab in Pigeon is just a directory under the root directory that contains
* a top-level `CMakeLists.txt`;
* a `pic.hpp`, where template parameters are specified that affect the builds of every module;
* a `pic_impl.hpp`, where details specific to this lab are provided.
Examples of these files are available under `examples/`. For most users, the `examples/CMakeLists.txt` and `examples/pic.hpp` can be directly copied into the lab for use with little to none modification. On the other hand, `pic_impl.hpp` is on a per-lab basis, so users must provide their own details. For instance, the physics setup to reproduce the earlier video can be found in `examples/pic_impl_pulsar.hpp`. For another, `examples/pic_impl_skeleton.hpp` provides the bare minimum to create an empty simulation that builds.

`pic_impl.hpp` should be the go-to file for all customizations. Customizations beyond that will have to involve modifying the underlying codebase, and therefore are considered advanced. Should such need occur, Pigeon supports shadowing a certain file in the codebase by a user provided alternative. Enable shadowing from `CMakeLists.txt`, make a directory named `shadows` (or whatever you chose to call in the `CMakeLists.txt`), and drop the shadowing file in it keeping the same filename and hierarchy. For example, to shadow `kernel/particle/migration_impl.hpp`, go to your label directory, run `mkdir -p shadows/particle`, then put your `migration_impl.hpp` there. (See, however, [It comes with a "secretary"'](#pgn) for automation)

Shadowing avoids directly modifying the codebase, which is important if you don't want the modification to affect all your project labs.

<a name="pgn"></a>
## It comes with a "secretary"! (Workflow automation)
Understandably, manually going over creating a lab directory, place essential files there takes many steps. However this is not even close to the end of the list. To greatly reduce this overhead, Pigeon shipsships with a very-high-level efficiency tool named `pgn`, a command-line tool designed
specifically for workflow on computer clusters as well as local laptops. (**NOTE** the following uses `pgn` to denote a path to the
*`your_Pigeon_root_dir`*`/pgn`
. Alternatively, create an alias in your shell rc file.)

The following steps are identified on a computer cluster if one starts from scratch and wants to obtain the pulsar result.
1. create a lab
2. build the simulation and if first time, run cmake before building
3. submit a job which runs the simulation to cluster scheduler
4. note down the jobID returned by the scheduler and if needed, a few words about what this run is about
5. when the run starts, note down the runtime data directory and associate it with the jobID
5. when the run finishes, look at the Pigeon generated report as well as `stdout` and `stderr`
6. resume the run as needed.

`pgn` streamlines these routines, which easily bypasses the tedium of maintaining
information integrity and coherency across various softwares. (Currently only PBS
and Slurm are supported.)
1. At Pigeon root directory, run `pgn new`*`lab_name`* to create a lab
2. `cd`*`lab_name`*
3. `pgn cb`, equivalent to first `pgn c`(runs CMake only) followed by `pgn b`(runs builds only). `pgn --clean` starts a clean build again. **Note** that you don't need to explictly create a build directory for `cmake`; the out-of-source build is realized by another project of mine called [Buick](https://github.com/hoorayphyer/Buick)
4. `pgn br -c `*`conf.toml`* to "run" the simulation, which is equivalent to `pgn b` followed by `pgn r`. The repeated `b` here is strongly recommended as it makes sure the most recent updates in your `pic_impl.hpp` get compiled and show up in the coming run. Actually quite a lot is happening behind `pgn r`. See [behind `pgn r`](#behind_pgn_r) below for a detailed explanation.
5. when run finishes, `pgn vw` to simultaneously open the report, the `jobID.OU` and the `jobID.ER` for quick review.
6. To resume a run, `pgn br --resume `*`path_to_checkpoint`*` -c `*`conf.toml`*.

<a name="behind_pgn_r"></a>
### behind `pgn r`
#### journaling : gathering the scattered
When launching simulations with `pgn r`, you will be prompted to enter a message about the run in `vim`. The message will be saved into the `journal.txt` file. The `journal.txt` is part of the journaling system of Pigeon, which aims exactly to keep information all in one place. For example, `pgn r` finishes by actually submitting a job, the jobID returned by the cluster job scheduler is also captured into `journal.txt`. Later when the simulation starts, the runtime data directory will be automatically saved to this file as well.
#### staging : the key to multi-tasking
Computer clusters usually allows submitting multiple jobs at one time. But caution must be taken if one wants to run multiple `pgn br -c `*`conf.toml`* inside one lab, because the cluster scheduler don't save copies of files pertaining to that run, such as `pic_impl.hpp` and your toml config if existed. That means, a launched run is potentially vulnerable to any changes you make in the lab. A reliable and carefree solution is to copy these essential files to a different path and submit the job there, with only one obvious drawback : tedious and error-prone if done by a human. Pigeon provides this automation through its staging system. The staging directory is also captured in `journal.txt`, as it's duty-bound to.
