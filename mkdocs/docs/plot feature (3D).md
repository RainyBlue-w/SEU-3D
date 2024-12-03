# Plot feature (3D)

<img src="https://pic-md-1259550128.cos.ap-nanjing.myqcloud.com/image-20241203103443025.png" alt="image-20241203103443025" style="zoom:50%;" />

​	In this tab the digital embryos in a 3D coordinates system are displayed, and they are formed by registering and aligning multiple layers of slices.

​	On the left is a accordion that contains some functions and configuration items.  On the middle are two 3D scatter with ***synchronized viewing angles***, which show gene expression and cell-type annotation respectively. On the right is a `cell-type selector`.

### Accordion on the left

- Select data

  Use  `Select data` to choose embryo in your interested stage, and the `germ layer` could also be specified to show part of the embryo.

<img src="https://pic-md-1259550128.cos.ap-nanjing.myqcloud.com/image-20241203104406456.png" alt="image-20241203104406456" style="zoom:50%;" />



- Plot options

  There are 3 subtabs in the `Plot options`, they are `Settings`, `Single` and `Multiple` respectively.

  - settings

    In `settings` you could configure the drawing options, which may be useful for generating the image you want.

  ​                                   <img src="https://pic-md-1259550128.cos.ap-nanjing.myqcloud.com/image-20241203104750979.png" alt="image-20241203104750979" style="zoom:50%;" /><img src="https://pic-md-1259550128.cos.ap-nanjing.myqcloud.com/image-20241203105322955.png" alt="image-20241203105322955" style="zoom:50%;" />

  - Single

    In `Single` you could choose your interested gene to plot on the 3D scatter. By clicking the color dot on the right, you could specify drawing color of the gene.

    ​                              <img src="https://pic-md-1259550128.cos.ap-nanjing.myqcloud.com/image-20241203105550325.png" alt="image-20241203105550325" style="zoom:50%;" /><img src="https://pic-md-1259550128.cos.ap-nanjing.myqcloud.com/image-20241203110027051.png" alt="image-20241203110027051" style="zoom:50%;" />

  - Multiple

    In `Multiple`, several genes can be displayed together on one 3D scatter. Color of each cell is assigned by mixing the specified color of selected gene, based on the expressing proportion in the cell.

<img src="C:\Users\慎司\AppData\Roaming\Typora\typora-user-images\image-20241203110813714.png" alt="image-20241203110813714" style="zoom:50%;" />

- Slicer

  ​	`Slicer` is used to subset part of the embryo. Check the `Preview` switch will show a black box  covering the range that cells in which will be preserved.

<img src="https://pic-md-1259550128.cos.ap-nanjing.myqcloud.com/image-20241203121435141.png" alt="image-20241203121435141" style="zoom:50%;" />



- Compute SVG (Moran)

  It's an analysis tool to find spatially variable genes in current cells ( filtered by `germ layer`, `Slicer` and `cell-type selector`). Results will be preserved in a drawer, and can be checked again by pressing `Result` button.

  <img src="https://pic-md-1259550128.cos.ap-nanjing.myqcloud.com/image-20241203121932141.png" alt="image-20241203121932141" style="zoom:50%;" />

