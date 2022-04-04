# My-rendering-engine
This is my software rendering engine using rasterisation, written in Rust

## Features
- 3D movement
- 3D voxel chunk collision, meshing and rendering with triangles
- textures can be animated
- tools to add colored lights, break and place blocks
- saving and loading chunk files
- day/night cycle
- rough chunk-based frustrum culling

## Problems
- a lot of things are hardcoded right now (default textures, font, window size...) so hard to customize
- there isn't rotation on every axis, and rotation of meshes can be buggy at times
- memory usage can be an issue
- not multi threaded, so still leaving a lot of performance on the table
