package ru.inspirit.fluid 
{

	import flash.utils.ByteArray;
	import com.joa_ebert.apparat.memory.Memory;
	import flash.geom.Matrix3D;
	import flash.display.BitmapData;
	/**
	 * @author Eugene Zatepyakin
	 */
	public class FluidRenderer 
	{
		public var solver:FluidSolver3D;
		
		public var verts:Vector.<Number>;
		
		public var pixies:FluidPix;
		
		public var posZBuffer:int;
		public var posScreen:int;
		public var screenBuffer:int = (780*460) << 2;
		public var buffer:ByteArray;
		
		public function FluidRenderer(solver:FluidSolver3D)
		{
			this.solver = solver;
			
			buffer = solver.bytes;
			
			var sizeX:int = solver.sizeX;
			var sizeY:int = solver.sizeY;
			var sizeZ:int = solver.sizeZ;
			
			var hsX:Number = sizeX * 0.5;
			var hsY:Number = sizeY * 0.5;
			var hsZ:Number = sizeZ * 0.5;
			
			var sizeX2:int = sizeX + 2;
			var sizeY2:int = sizeY + 2;
			
			var i:int, j:int, k:int;
			var ind:int;
			var scale:Number = 4.5;//780 / 4 / sizeX;
			
			var p:FluidPix;
			
			for(k = sizeZ; k > 0; --k) 
			{
		        for(j = sizeY; j > 0; --j) 
		        {
		        	ind = ((k * sizeX2 * sizeY2) + (j * sizeX2) + sizeX) << 2;
		            for(i = sizeX; i > 0; --i) 
		            {		            	
		            	if( null == p )
		            	{
							p = pixies = new FluidPix();
		            	} else {
							p = p.next = new FluidPix();
						}
						p.x = (i - hsX) * scale;
						p.y = (j - hsY) * scale;
						p.z = (k - hsZ) * scale;
						
						p.ind = ind;
						p.ir = solver.posR + ind;
						p.ig = solver.posG + ind;
						p.ib = solver.posB + ind;
		            	
		            	ind -= 4;
		            }
		        }
			}
			
			posZBuffer = solver.posZBuffer;
			posScreen = solver.posScreen;
		}

		public function render(bmp:BitmapData, matrix:Matrix3D, focalLength:Number = 60):void
		{
			var pz:Number;
			
			var p00: Number = matrix.rawData[ 0x0 ];
			var p01: Number = matrix.rawData[ 0x1 ];
			var p02: Number = matrix.rawData[ 0x2 ];
			var p10: Number = matrix.rawData[ 0x4 ];
			var p11: Number = matrix.rawData[ 0x5 ];
			var p12: Number = matrix.rawData[ 0x6 ];
			var p20: Number = matrix.rawData[ 0x8 ];
			var p21: Number = matrix.rawData[ 0x9 ];
			var p22: Number = matrix.rawData[ 0xa ];
			var p32: Number = matrix.rawData[ 0xe ];
			
			var bw:int = bmp.width;
			var bh:int = bmp.height;
			var cx: Number = bw * 0.5;
			var cy: Number = bh * 0.5;
			var w:Number, x:Number, y:Number, z:Number;
			var i:int, j:int, ind:int;
			var c:int;
			
			var p:FluidPix = pixies;
			i = posZBuffer;
			j = posZBuffer + screenBuffer;
			while(i < j) 
			{
				Memory.writeFloat(-0xFFFFFF, i);
				Memory.writeInt(0, i+screenBuffer);
				i += 4;
			}
			
			var screenOffset:int = posScreen;
			var bufferWidth:int = bw << 2;
			var bufferMax:int = ((bw * bh) << 2) - 4;
			var bufferMin:int = -1;
			var bufferIndex:int;
			
			while(p)
			{
				x = p.x; y = p.y; z = p.z;
				pz = focalLength + x * p02 + y * p12 + z * p22 + p32;
				c = int(Memory.readFloat(p.ir) * 0xFF) << 16 | int(Memory.readFloat(p.ig) * 0xFF) << 8 | int(Memory.readFloat(p.ib) * 0xFF);
				
				if(pz > 0 && c > 0)
				{
					i = int( ( w = focalLength / pz ) * ( x * p00 + y * p10 + z * p20 ) + cx );
					j = int( w * ( x * p01 + y * p11 + z * p21 ) + cy );
					
					if(bufferMin < ( bufferIndex = int( ( i << 2 ) + int( j * bufferWidth ) ) ) && bufferIndex < bufferMax
						&& !(w < Memory.readFloat( (ind = posZBuffer + bufferIndex) )))
					{
						//var prevC:int = Memory.readInt( screenOffset + bufferIndex );
						/*c = (0xFF0000 & ((c & 0xFF0000) +
										      (((int(prevC & 0xFF0000) -
										      int(c & 0xFF0000)) * 0x00) >>8))) |
										   (0x00FF00 & ((c & 0x00FF00) + (((int(prevC & 0x00FF00) -
										      int(c & 0x00FF00)) * 0x00) >>8))) |
										   (0x0000FF & ((c & 0x0000FF) +
										      (((int(prevC & 0x0000FF) -
										      int(c & 0x0000FF)) * 0x00) >>8)));*/

						Memory.writeInt( c, screenOffset + bufferIndex );
						Memory.writeFloat(w, ind);
					}
				}
				
				p = p.next;
			}
			buffer.position = screenOffset;
			bmp.setPixels(bmp.rect, buffer);
		}
	}
}

internal final class FluidPix
{
	public var next:FluidPix;
	
	public var x:Number;
	public var y:Number;
	public var z:Number;
	
	public var ind:int;
	
	public var ir:int;
	public var ig:int;
	public var ib:int;
}
