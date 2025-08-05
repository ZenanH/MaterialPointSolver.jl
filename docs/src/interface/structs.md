# Structs design

!!! note

    Here is an interactive graph to illustrate the basic structs in this solver.
    
```@raw html
<iframe frameborder="0" style="width:100%;height:993px;" src="https://viewer.diagrams.net/?lightbox=1&target=blank&highlight=0000ff&nav=1&title=types.drawio&transparent=1&dark=1#Uhttps%3A%2F%2Fdrive.google.com%2Fuc%3Fid%3D12XsHGS7FoGt8cuzTLdUagSsZHNrODts8%26export%3Ddownload" allowtransparency="true"></iframe>
```

Green blocks are structs that users need to interact with, and they follow a similar hierarchical logic. To maximize support for users to easily implement custom algorithms, we have added an Extra structure to each structure. For example, for `Particle`, it is `UserParticleExtra`.

!!! note

    Although `Args` has a user-customizable struct called `UserArgsExtra`, like other `Args`, it cannot be directly used as input for the kernel function; it can only pass its values (bits) to the kernel function.

When the user instantiates green blocks, the data will be automatically transferred between different computing backends based on the user's selection. For internal information about the types, please refer to the source code.