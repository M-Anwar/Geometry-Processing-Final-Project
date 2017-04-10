# Implicit Skinning: Deformation with Cloth and Skin Contact Modelling
## CSC2521: Geometry Processing Final Project

> **To get started:** Fork this repository then issue
> 
>     git clone --recursive http://github.com/[username]/Geometry-Processing-Final-Project.git
>

## Installation and Compilation

See
[introduction](http://github.com/alecjacobson/geometry-processing-introduction).

## Execution

Once built, you can execute the assignment from inside the `build/` by running:

    ./Implicit_Skinning

### Implicit Skinning Example
When running the implicit skinning example, the following key commands will modify viewer to show different visualizations.

* 'w' -> show skinning weights
* 'e' -> switch bone index for different bone visualizations
* 'n' -> show normals for entire mesh
* 'p' -> Show an individual partitions, press 'e' to cycle between bones
* 'm' -> Show the sampled HRBF points used to approximate the surface
* 'b' -> Show individual HRBF functions on each bone
* 'v' -> Show ALL HRBF functions for entire mesh.
* ' '(space) -> Begin animation
* 'd' -> Switch between dual quaternion and implicit skinning.

### Cloth Simulation Example
When running the cloth simulation examplle, the following key commands will change the behaviour:
* 'd' -> Enable time step (run the simulation)
* ' '(space)-> Enable animation (automatically run the time step).
* 'a' -> Manual time step **BUGGY NOT IMPLEMENTED FULLY**

### Switching Examples
To switch between examples go to the main.cpp file and locate the following variable:

```C++
  bool run_cloth_sim = true;
```
setting to **true** will run the cloth simulation example, and **false** will run the implicit skinning example.

Yes, that is a bit cumbersome, but if you just want to see the final result you can look below for some GIFS and figures. The cloth simulation is still a bit wonky and will fail for very fast moving objects, so you may need to adjust the time-step.


## Final Results

### Implicit Skinning Example
The following shows the final result between implicit skinning and dual quaternion. The GIFS are not syncronized because my implementation is quite a bit slower then DQS but the results are as expected. Implicit skinning provides a more natural looking deformation for small to medium bending angles. It avoids the unnatural 'bulges' that appear in the dual quaternion solution. However, at extreme bending angles, both methods result in significant self intersections. 

Dual Quaternion Skinning   |  Implicit Skinning
:-------------------------:|:-------------------------:
![](images/dual_quat.gif)  |  ![](images/implicit.gif)

### Cloth Simulation Example

The following are the final results from the cloth simulation collision example. Here we use the computed implicit surface representation of the underlying mesh to perform efficient collision detection against the cloth. The results are not perfect, the hand pokes through the cloth because in this example the hands did not have implicit surfaces associated with them. Furthermore, sometimes the underlying mesh pokes through the cloth and the collisions are not resolved correctly. The overall solution is not that great and some more refinement is definitely required.

Cloth Simulation with No Collision   |  Cloth Simulation with Implicit Collisions
:-------------------------:|:-------------------------:
![](images/cloth_sim_nocol.png)  |  ![](images/cloth_sim_final.png)

Live result of cloth collision detection. The red dots indicate collision points evaluated using the implicit surface distance functions. In this example only the forarm is being used to generate the collision points.  

![](images/cloth_sim_collision.gif)


