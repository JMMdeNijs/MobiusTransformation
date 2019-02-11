# MobiusTransformation
This  C++ generates a Mobius transformation of an image. 

An elegant explanation of the Mobius transformation can be found: https://www.youtube.com/watch?v=0z1fIsUNhO4

The program comprises two functional parts:
a)	a basic Mobius Transformation to an new image of double width and height
b)	filling up the missing pixels in the new image

code explores different methods to fill up the missing pixels. Using key boards commands, the method can be changed:

	'B': no filling but the bare transformation
	'Q': the quads spanned by neighboring are split in 2 or 4 subquads and the RGB color of the new corners are calculated and added in the image. This proces works recursively.
	'T': the quads spanned by neighboring are split in 2 triangles and the RGB color of the inside pixels are calculated and added in the image,
	'P': Pause of the rotation

	Up and Down arrows can be used to fill neigbor pixels with RGB = {0, 0, 0}. 
	
	
An illustration can be found: https://www.dropbox.com/l/scl/AAAmcJUaM5Ly0HgrzmExK6r3J35uyvPpFZo


