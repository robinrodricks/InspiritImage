package ru.inspirit.fluid 
{

	import com.joa_ebert.apparat.memory.Memory;
	import flash.utils.Endian;
	import flash.utils.ByteArray;
	/**
	 * @author Eugene Zatepyakin
	 */
	public final class FluidSolver3D 
	{
		public var sizeX:int;
		public var sizeY:int;
		public var sizeZ:int;
		
		protected var sizeX2:int;
		protected var sizeY2:int;
		protected var sizeZ2:int;
		protected var sizeXY2:int;
		
		public var dt:Number = 0.5;
		public var colorDiffusion:Number = 0.0001;
		public var viscosity:Number = 0.0001;
		public var fadeSpeed:Number = 0.007;
		public var solverIter:int = 4;
		
		public const bytes:ByteArray = new ByteArray();
		protected var bytesLength:int;
		protected var bytesPerArray:int;
		
		public var posR:int;
		public var posG:int;
		public var posB:int;
		public var posPrevR:int;
		public var posPrevG:int;
		public var posPrevB:int;
		
		public var posVelX:int;
		public var posVelY:int;
		public var posVelZ:int;
		public var posPrevVelX:int;
		public var posPrevVelY:int;
		public var posPrevVelZ:int;
		public var posZBuffer:int;
		public var posScreen:int;
		
		public var binIndex:Vector.<FluidBIN>;
		public var bins:FluidBIN;
		public var binsFull:FluidBIN;
		public var binsYX:FluidBIN;
		public var binsZX:FluidBIN;
		public var binsZY:FluidBIN;
		
		public function FluidSolver3D(width:int = 20, height:int = 20, depth:int = 20)
		{
			sizeX = width;
			sizeY = height;
			sizeZ = depth;
			
			sizeX2 = width + 2;
			sizeY2 = height + 2;
			sizeZ2 = depth + 2;
			
			sizeXY2 = sizeX2 * sizeY2;
			
			reset();
			buildBins();
		}
		
		public function reset():void
		{
			bytesPerArray = (sizeX2 * sizeY2 * sizeZ2) << 2;
			
			posR = 0;
			posG = bytesPerArray;
			posB = bytesPerArray * 2;
			posPrevR = bytesPerArray * 3;
			posPrevG = bytesPerArray * 4;
			posPrevB = bytesPerArray * 5;
			
			posVelX = bytesPerArray * 6;
			posVelY = bytesPerArray * 7;
			posVelZ = bytesPerArray * 8;
			posPrevVelX = bytesPerArray * 9;
			posPrevVelY = bytesPerArray * 10;
			posPrevVelZ = bytesPerArray * 11;
			posZBuffer = bytesPerArray * 12;
			posScreen = posZBuffer + ((780*460) << 2);
			
			bytesLength = posScreen + ((780*460) << 2);
			
			bytes.length = 0;
			bytes.endian = Endian.LITTLE_ENDIAN;
			bytes.length = bytesLength;
			
			Memory.select(bytes);
		}
		
		public function update():void
		{
			var bin:FluidBIN;			
			
			// diffuse XYZ			
			var a:Number = dt * viscosity * sizeX * sizeY;
			solveLinearXYZ(a, 1 / (a * 6 + 1));
			
			project1(posPrevVelX, posPrevVelY, posPrevVelZ, posVelX, posVelY);
			
			advectXYZ(dt);
			
			project2(posVelX, posVelY, posVelZ, posPrevVelX, posPrevVelY);
			
			if( colorDiffusion != 0 && dt != 0 )
            {
				a = dt * colorDiffusion * sizeX * sizeY;
            	solveLinearRGB(a, 1 / (a * 6 + 1));
            } else
            {
            	bin = bins;
            	while(bin)
            	{
            		Memory.writeFloat(Memory.readFloat(bin.ri), bin.pri);
            		Memory.writeFloat(Memory.readFloat(bin.gi), bin.pgi);
            		Memory.writeFloat(Memory.readFloat(bin.bi), bin.pbi);
            		bin = bin.next;
				}
			}
            
            advectRGB(dt);
            
            // Fade RGB
            var holdAmount:Number = 1 - fadeSpeed;
			bin = binsFull;
			while(bin)
        	{
        		Memory.writeFloat(Memory.readFloat(bin.ri) * holdAmount, bin.ri);
        		Memory.writeFloat(Memory.readFloat(bin.gi) * holdAmount, bin.gi);
        		Memory.writeFloat(Memory.readFloat(bin.bi) * holdAmount, bin.bi);
        		bin = bin.next;
			}
		}
		
		protected function advectXYZ(dt:Number):void
		{
			var i0:int, i1:int, j0:int, j1:int, k0:int, k1:int;

			var dtx:Number = dt * sizeX;
			var dty:Number = dt * sizeY;
			var dtz:Number = dt * sizeZ;
			
			var s0:Number, s1:Number, t0:Number, t1:Number, u0:Number, u1:Number;
			var x:Number, y:Number, z:Number;
			
			var ind2:int;
			
			var zdem:int = sizeXY2 << 2;
			var xdem:int = sizeX2 << 2;
			
			var bin:FluidBIN = bins;

			while(bin)
			{
						x = bin.i - dtx * Memory.readFloat(bin.pxi);
		                y = bin.j - dty * Memory.readFloat(bin.pyi);
		                z = bin.k - dtz * Memory.readFloat(bin.pzi);
		
		                if(x < 0.5) x = 0.5;
						if(x > sizeX + 0.5) x = sizeX + 0.5;
		                i0 = (x);
		                i1 = i0 + 1;
		                if(y < 0.5) y = 0.5;
		                if(y > sizeY + 0.5) y = sizeY + 0.5;
		                j0 = (y);
		                j1 = j0 + 1;
		                if(z < 0.5) z = 0.5;
		                if(z > sizeZ + 0.5) z = sizeZ + 0.5;
		                k0 = (z);
		                k1 = k0 + 1;
		
		                s1 = x - i0;
		                s0 = 1.0 - s1;
		                t1 = y - j0;
		                t0 = 1.0 - t1;
		                u1 = z - k0;
		                u0 = 1.0 - u1;
						
		                ind2 = ((k0 * sizeXY2) + (j0 * sizeX2) + i0) << 2;
		                
		                var ii0:int = ind2 + zdem;
		                var ii1:int = ind2 + xdem;
		                var ii2:int = ii1 + zdem;
		                var ii3:int = ind2 + 4;
		                var ii4:int = ii3 + zdem;
		                var ii5:int = ind2 + xdem + 4;
		                var ii6:int = ind2 + xdem + 4 + zdem;
		                
		                Memory.writeFloat(
		                    s0 * ( t0 * (u0 * Memory.readFloat(posPrevVelX + ind2)
		                                +u1 * Memory.readFloat(posPrevVelX + ii0))
		                        +( t1 * (u0 * Memory.readFloat(posPrevVelX + ii1)
		                                +u1 * Memory.readFloat(posPrevVelX + ii2))))
		                   +s1 * ( t0 * (u0 * Memory.readFloat(posPrevVelX + ii3)
		                                +u1 * Memory.readFloat(posPrevVelX + ii4))
		                        +( t1 * (u0 * Memory.readFloat(posPrevVelX + ii5)
		                                +u1 * Memory.readFloat(posPrevVelX + ii6)))), 
		                                bin.xi);
						Memory.writeFloat(
		                    s0 * ( t0 * (u0 * Memory.readFloat(posPrevVelY + ind2)
		                                +u1 * Memory.readFloat(posPrevVelY + ii0))
		                        +( t1 * (u0 * Memory.readFloat(posPrevVelY + ii1)
		                                +u1 * Memory.readFloat(posPrevVelY + ii2))))
		                   +s1 * ( t0 * (u0 * Memory.readFloat(posPrevVelY + ii3)
		                                +u1 * Memory.readFloat(posPrevVelY + ii4))
		                        +( t1 * (u0 * Memory.readFloat(posPrevVelY + ii5)
		                                +u1 * Memory.readFloat(posPrevVelY + ii6)))), 
		                                bin.yi);
						Memory.writeFloat(
		                   s0 * ( t0 * (u0 * Memory.readFloat(posPrevVelZ + ind2)
		                                +u1 * Memory.readFloat(posPrevVelZ + ii0))
		                        +( t1 * (u0 * Memory.readFloat(posPrevVelZ + ii1)
		                                +u1 * Memory.readFloat(posPrevVelZ + ii2))))
		                   +s1 * ( t0 * (u0 * Memory.readFloat(posPrevVelZ + ii3)
		                                +u1 * Memory.readFloat(posPrevVelZ + ii4))
		                        +( t1 * (u0 * Memory.readFloat(posPrevVelZ + ii5)
		                                +u1 * Memory.readFloat(posPrevVelZ + ii6)))), 
		                                bin.zi);
				bin = bin.next;
			}
			
			//keepInBounds(1, posVelX);
			//keepInBounds(2, posVelY);
			//keepInBounds(3, posVelZ);
			keepInBoundsXYZ(posVelX, posVelY, posVelZ);
		}
		
		protected function advectRGB(dt:Number):void
		{
			var i0:int, i1:int, j0:int, j1:int, k0:int, k1:int;

			var dtx:Number = dt * sizeX;
			var dty:Number = dt * sizeY;
			var dtz:Number = dt * sizeZ;
			
			var s0:Number, s1:Number, t0:Number, t1:Number, u0:Number, u1:Number;
			var x:Number, y:Number, z:Number;
			
			var ind2:int;
			
			var zdem:int = sizeXY2 << 2;
			var xdem:int = sizeX2 << 2;
			
			var bin:FluidBIN = bins;

			while(bin)
			{
						x = bin.i - dtx * Memory.readFloat(bin.xi);
		                y = bin.j - dty * Memory.readFloat(bin.yi);
		                z = bin.k - dtz * Memory.readFloat(bin.zi);
		
		                if(x < 0.5) x = 0.5;
						if(x > sizeX + 0.5) x = sizeX + 0.5;
		                i0 = (x);
		                i1 = i0 + 1;
		                if(y < 0.5) y = 0.5;
		                if(y > sizeY + 0.5) y = sizeY + 0.5;
		                j0 = (y);
		                j1 = j0 + 1;
		                if(z < 0.5) z = 0.5;
		                if(z > sizeZ + 0.5) z = sizeZ + 0.5;
		                k0 = (z);
		                k1 = k0 + 1;
		
		                s1 = x - i0;
		                s0 = 1.0 - s1;
		                t1 = y - j0;
		                t0 = 1.0 - t1;
		                u1 = z - k0;
		                u0 = 1.0 - u1;
						
		                ind2 = ((k0 * sizeXY2) + (j0 * sizeX2) + i0) << 2;
		                
		                var ii0:int = ind2 + zdem;
		                var ii1:int = ind2 + xdem;
		                var ii2:int = ii1 + zdem;
		                var ii3:int = ind2 + 4;
		                var ii4:int = ii3 + zdem;
		                var ii5:int = ind2 + xdem + 4;
		                var ii6:int = ind2 + xdem + 4 + zdem;
		                
		                Memory.writeFloat(
		                    s0 * ( t0 * (u0 * Memory.readFloat(posPrevR + ind2)
		                                +u1 * Memory.readFloat(posPrevR + ii0))
		                        +( t1 * (u0 * Memory.readFloat(posPrevR + ii1)
		                                +u1 * Memory.readFloat(posPrevR + ii2))))
		                   +s1 * ( t0 * (u0 * Memory.readFloat(posPrevR + ii3)
		                                +u1 * Memory.readFloat(posPrevR + ii4))
		                        +( t1 * (u0 * Memory.readFloat(posPrevR + ii5)
		                                +u1 * Memory.readFloat(posPrevR + ii6)))), 
		                                bin.ri);
						Memory.writeFloat(
		                    s0 * ( t0 * (u0 * Memory.readFloat(posPrevG + ind2)
		                                +u1 * Memory.readFloat(posPrevG + ii0))
		                        +( t1 * (u0 * Memory.readFloat(posPrevG + ii1)
		                                +u1 * Memory.readFloat(posPrevG + ii2))))
		                   +s1 * ( t0 * (u0 * Memory.readFloat(posPrevG + ii3)
		                                +u1 * Memory.readFloat(posPrevG + ii4))
		                        +( t1 * (u0 * Memory.readFloat(posPrevG + ii5)
		                                +u1 * Memory.readFloat(posPrevG + ii6)))), 
		                                bin.gi);
						Memory.writeFloat(
		                   s0 * ( t0 * (u0 * Memory.readFloat(posPrevB + ind2)
		                                +u1 * Memory.readFloat(posPrevB + ii0))
		                        +( t1 * (u0 * Memory.readFloat(posPrevB + ii1)
		                                +u1 * Memory.readFloat(posPrevB + ii2))))
		                   +s1 * ( t0 * (u0 * Memory.readFloat(posPrevB + ii3)
		                                +u1 * Memory.readFloat(posPrevB + ii4))
		                        +( t1 * (u0 * Memory.readFloat(posPrevB + ii5)
		                                +u1 * Memory.readFloat(posPrevB + ii6)))), 
		                                bin.bi);
		    	bin = bin.next;                           
			}
		    
		    //keepInBounds(0, posR);
			//keepInBounds(0, posG);
			//keepInBounds(0, posB);
		}
		
		protected function project2(velX:int, velY:int, velZ:int, pPos:int, divPos:int):void
		{
			var xdem:int = sizeX2 << 2;
			var zdem:int = sizeXY2 << 2;
			var h:Number = -0.5 / sizeX;
			
			var bin:FluidBIN = bins;

			while(bin)
			{
				Memory.writeFloat(
									h *
									(Memory.readFloat(bin.xi+4)
									-Memory.readFloat(bin.xi-4)
									+Memory.readFloat(bin.yi+xdem)
									-Memory.readFloat(bin.yi-xdem)
									+Memory.readFloat(bin.zi+zdem)
									-Memory.readFloat(bin.zi-zdem)
									), 
									bin.pyi);
				Memory.writeFloat(0, bin.pxi);
				
				bin = bin.next;
			}
			
			keepInBounds(0, divPos);
			keepInBounds(0, pPos);
			
			var c:Number = 1/6;
			var k:int, ind:int;
			for (k = 0; k < solverIter; ++k) 
			{
				bin = bins;
				while(bin)
		    	{
		    		ind = bin.pxi;
						Memory.writeFloat(
		            						(Memory.readFloat(bin.pyi)
		            						+ (Memory.readFloat(ind + 4)
		            								+Memory.readFloat(ind - 4)
		            								+Memory.readFloat(ind + xdem)
		            								+Memory.readFloat(ind - xdem)
		            								+Memory.readFloat(ind + zdem)
		            								+Memory.readFloat(ind - zdem)
		            								)) * c
		            								, ind);
		            	bin = bin.next;
		    	}
				keepInBounds(0, pPos);
			}
			
		    bin = bins;
		    h = 0.5 * sizeX;
		    while(bin)
		    {
		    	ind = bin.pxi;
		            	Memory.writeFloat(
		            						Memory.readFloat(bin.xi)
		            						-(h * 
		            						(Memory.readFloat(ind+4)
		            						-Memory.readFloat(ind-4))), 
		            						bin.xi);
						Memory.writeFloat(
		            						Memory.readFloat(bin.yi)
		            						-(h * 
		            						(Memory.readFloat(ind + xdem)
		            						-Memory.readFloat(ind - xdem))), 
		            						bin.yi);
		            	Memory.writeFloat(
		            						Memory.readFloat(bin.zi)
		            						-(h * 
		            						(Memory.readFloat(ind + zdem)
		            						-Memory.readFloat(ind - zdem))), 
		            						bin.zi);
		    	bin = bin.next;
		    }
		    
		    //keepInBounds(1, velX);
			//keepInBounds(2, velY);
			//keepInBounds(3, velZ);
			keepInBoundsXYZ(velX, velY, velZ);
		}
		
		protected function project1(velX:int, velY:int, velZ:int, pPos:int, divPos:int):void
		{
			var xdem:int = sizeX2 << 2;
			var zdem:int = sizeXY2 << 2;
			var h:Number = -0.5 / sizeX;
			
			var bin:FluidBIN = bins;

			while(bin)
			{
				Memory.writeFloat(
									h *
									(Memory.readFloat(bin.pxi+4)
									-Memory.readFloat(bin.pxi-4)
									+Memory.readFloat(bin.pyi+xdem)
									-Memory.readFloat(bin.pyi-xdem)
									+Memory.readFloat(bin.pzi+zdem)
									-Memory.readFloat(bin.pzi-zdem)
									), 
									bin.yi);
				Memory.writeFloat(0, bin.xi);
				
				bin = bin.next;
			}
			
			keepInBounds(0, divPos);
			keepInBounds(0, pPos);
			
			var c:Number = 1/6;
			var k:int, ind:int;
			for (k = 0; k < solverIter; ++k) 
			{
				bin = bins;
				while(bin)
		    	{
		    		ind = bin.xi;
						Memory.writeFloat(
		            						(Memory.readFloat(bin.yi)
		            						+ (Memory.readFloat(ind + 4)
		            								+Memory.readFloat(ind - 4)
		            								+Memory.readFloat(ind + xdem)
		            								+Memory.readFloat(ind - xdem)
		            								+Memory.readFloat(ind + zdem)
		            								+Memory.readFloat(ind - zdem)
		            								)) * c
		            								, ind);
		            	bin = bin.next;
		    	}
				keepInBounds(0, pPos);
			}
			
		    bin = bins;
		    h = 0.5 * sizeX;
		    while(bin)
		    {
		    	ind = bin.xi;
		            	Memory.writeFloat(
		            						Memory.readFloat(bin.pxi)
		            						-(h * 
		            						(Memory.readFloat(ind+4)
		            						-Memory.readFloat(ind-4))), 
		            						bin.pxi);
						Memory.writeFloat(
		            						Memory.readFloat(bin.pyi)
		            						-(h * 
		            						(Memory.readFloat(ind + xdem)
		            						-Memory.readFloat(ind - xdem))), 
		            						bin.pyi);
		            	Memory.writeFloat(
		            						Memory.readFloat(bin.pzi)
		            						-(h * 
		            						(Memory.readFloat(ind + zdem)
		            						-Memory.readFloat(ind - zdem))), 
		            						bin.pzi);
		    	bin = bin.next;
		    }
		    
		    //keepInBounds(1, velX);
			//keepInBounds(2, velY);
			//keepInBounds(3, velZ);
			keepInBoundsXYZ(velX, velY, velZ);
		}
		
		protected function solveLinearXYZ(a:Number, c:Number):void
		{
			var k:int, ind:int;
			
			var xdem:int = sizeX2 << 2;
			var zdem:int = sizeXY2 << 2;
			
			var bin:FluidBIN;
			
			for (k = 0; k < solverIter; ++k) 
			{
				bin = bins;
				while(bin)
				{
					ind = bin.pxi;
	            	Memory.writeFloat(
	            						(Memory.readFloat(bin.xi)
	            						+ a * (Memory.readFloat(ind + 4)
	            								+Memory.readFloat(ind - 4)
	            								+Memory.readFloat(ind + xdem)
	            								+Memory.readFloat(ind - xdem)
	            								+Memory.readFloat(ind + zdem)
	            								+Memory.readFloat(ind - zdem)
	            								)) * c, 
	            								ind);
					ind = bin.pyi;
					Memory.writeFloat(
	            						(Memory.readFloat(bin.yi)
	            						+ a * (Memory.readFloat(ind + 4)
	            								+Memory.readFloat(ind - 4)
	            								+Memory.readFloat(ind + xdem)
	            								+Memory.readFloat(ind - xdem)
	            								+Memory.readFloat(ind + zdem)
	            								+Memory.readFloat(ind - zdem)
	            								)) * c, 
	            								ind);
					ind = bin.pzi;
	            	Memory.writeFloat(
	            						(Memory.readFloat(bin.zi)
	            						+ a * (Memory.readFloat(ind + 4)
	            								+Memory.readFloat(ind - 4)
	            								+Memory.readFloat(ind + xdem)
	            								+Memory.readFloat(ind - xdem)
	            								+Memory.readFloat(ind + zdem)
	            								+Memory.readFloat(ind - zdem)
	            								)) * c, 
	            								ind);
			        	bin = bin.next;
				}
				//keepInBounds(1, posPrevVelX);
				//keepInBounds(2, posPrevVelY);
				//keepInBounds(3, posPrevVelZ);
				keepInBoundsXYZ(posPrevVelX, posPrevVelY, posPrevVelZ);
			}
		}
		
		protected function solveLinearRGB(a:Number, c:Number):void
		{
			var k:int, ind:int;
			
			var xdem:int = sizeX2 << 2;
			var zdem:int = sizeXY2 << 2;
			
			var bin:FluidBIN;
			
			for (k = 0; k < solverIter; ++k) 
			{
				bin = bins;
				while(bin)
				{
			        ind = bin.pri;
	            	Memory.writeFloat(
	            						(Memory.readFloat(bin.ri)
	            						+ a * (Memory.readFloat(ind + 4)
	            								+Memory.readFloat(ind - 4)
	            								+Memory.readFloat(ind + xdem)
	            								+Memory.readFloat(ind - xdem)
	            								+Memory.readFloat(ind + zdem)
	            								+Memory.readFloat(ind - zdem)
	            								)) * c, 
	            								ind);
					ind = bin.pgi;
					Memory.writeFloat(
	            						(Memory.readFloat(bin.gi)
	            						+ a * (Memory.readFloat(ind + 4)
	            								+Memory.readFloat(ind - 4)
	            								+Memory.readFloat(ind + xdem)
	            								+Memory.readFloat(ind - xdem)
	            								+Memory.readFloat(ind + zdem)
	            								+Memory.readFloat(ind - zdem)
	            								)) * c, 
	            								ind);
					ind = bin.pbi;
	            	Memory.writeFloat(
	            						(Memory.readFloat(bin.bi)
	            						+ a * (Memory.readFloat(ind + 4)
	            								+Memory.readFloat(ind - 4)
	            								+Memory.readFloat(ind + xdem)
	            								+Memory.readFloat(ind - xdem)
	            								+Memory.readFloat(ind + zdem)
	            								+Memory.readFloat(ind - zdem)
	            								)) * c, 
	            								ind);
			        	bin = bin.next;
				}
				//keepInBounds(0, posPrevR);
				//keepInBounds(0, posPrevG);
				//keepInBounds(0, posPrevB);
			}
		}
		
		protected function keepInBounds(b:int, pos:int):void
		{
			var xdem:int = sizeX2 << 2;
			var xend:int = sizeX << 2;
			var yend:int = xdem * sizeY;
			var zdem:int = sizeXY2 << 2;
			var zend:int = zdem * sizeZ;
			var poszendpzdem:int = pos + zend + zdem;
			var posyendpxdem:int = pos + yend + xdem;
			var poszdem:int = pos + zdem;
			var poszend:int = pos + zend;
			var posxdem:int = pos + xdem;
			var posxendp4:int = pos + xend + 4;
			var posyend:int = pos + yend;
			var posxend:int = pos + xend;
			var pos4:int = pos + 4;
			
			var bin:FluidBIN;
			
			if(b == 3)
			{
				bin = binsYX;
				while(bin)
				{
					Memory.writeFloat(
				    						-Memory.readFloat(bin.ind + poszdem),
				    						bin.ind + pos);
				    Memory.writeFloat(
				    						-Memory.readFloat(bin.ind + poszend),
				    						bin.ind + poszendpzdem);
					bin = bin.nextYX;        
				}
			} 
			else 
			{
				bin = binsYX;
				while(bin)
				{
					Memory.writeFloat(
				    						Memory.readFloat(bin.ind + poszdem),
				    						bin.ind + pos);
					Memory.writeFloat(
				    						Memory.readFloat(bin.ind + poszend),
				    						bin.ind + poszendpzdem);
				    bin = bin.nextYX;
				}
			}
			
			if(b == 2)
			{
				bin = binsZX;
				while(bin)
				{
					Memory.writeFloat(
				    						-Memory.readFloat(bin.ind + posxdem),
				    						bin.ind + pos);
				    Memory.writeFloat(
				    						-Memory.readFloat(bin.ind + posyend),
				    						bin.ind + posyendpxdem);
					bin = bin.nextZX;        
				}
			}
			else
			{
				bin = binsZX;
				while(bin)
				{
					Memory.writeFloat(
				    						Memory.readFloat(bin.ind + posxdem),
				    						bin.ind + pos);
				    Memory.writeFloat(
				    						Memory.readFloat(bin.ind + posyend),
				    						bin.ind + posyendpxdem);
					bin = bin.nextZX;        
				}
			}
			
			if(b == 1)
			{
				bin = binsZY;
				while(bin)
				{
					Memory.writeFloat(
				    						-Memory.readFloat(bin.ind + pos4),
				    						bin.ind + pos);
				    Memory.writeFloat(
				    						-Memory.readFloat(bin.ind + posxend),
				    						bin.ind + posxendp4);
					bin = bin.nextZY;        
				}
			}
			else
			{
				bin = binsZY;
				while(bin)
				{
					Memory.writeFloat(
				    						Memory.readFloat(bin.ind + pos4),
				    						bin.ind + pos);
				    Memory.writeFloat(
				    						Memory.readFloat(bin.ind + posxend),
				    						bin.ind + posxendp4);
					bin = bin.nextZY;        
				}
			}
			
			Memory.writeFloat( 
								0.33 * (Memory.readFloat(pos4)
                                  + Memory.readFloat(posxdem)
                                  + Memory.readFloat(poszdem)),
                                  pos);
			Memory.writeFloat( 
								0.33 * (Memory.readFloat(pos+corner21)
                                  + Memory.readFloat(posyend)
                                  + Memory.readFloat(pos + corner23)),
                                  pos + corner24);
			Memory.writeFloat( 
								0.33 * (Memory.readFloat(pos+corner31)
                                  + Memory.readFloat(pos + corner32)
                                  + Memory.readFloat(poszend)),//
                                  pos + corner34);
			Memory.writeFloat( 
								0.33 * (Memory.readFloat(pos+corner41)//
                                  + Memory.readFloat(pos + corner42)
                                  + Memory.readFloat(pos + corner43)),
                                  pos + corner44);
			Memory.writeFloat( 
								0.33 * (Memory.readFloat(pos+xend)
                                  + Memory.readFloat(pos + corner52)
                                  + Memory.readFloat(pos + corner53)),
                                  pos + corner54);
			Memory.writeFloat( 
								0.33 * (Memory.readFloat(pos+corner61)
                                  + Memory.readFloat(pos + corner62)
                                  + Memory.readFloat(pos + corner63)),
                                  pos + corner64);
			Memory.writeFloat( 
								0.33 * (Memory.readFloat(pos+corner71)
                                  + Memory.readFloat(pos + corner72)
                                  + Memory.readFloat(pos + corner73)),
                                  pos + corner74);
			Memory.writeFloat( 
								0.33 * (Memory.readFloat(pos+corner81)
                                  + Memory.readFloat(pos + corner82)
                                  + Memory.readFloat(pos + corner83)),
                                  pos + corner84);
		}
		
		protected function keepInBoundsXYZ(posX:int, posY:int, posZ:int):void
		{
			var xdem:int = sizeX2 << 2;
			var xend:int = sizeX << 2;
			var yend:int = xdem * sizeY;
			var zdem:int = sizeXY2 << 2;
			var zend:int = zdem * sizeZ;
			var ind:int;
			
			var posXzdem:int = posX + zdem;
			var posXzend:int = posX + zend;
			var posXzendzdem:int = posX + zend + zdem;
			
			var posYzdem:int = posY + zdem;
			var posYzend:int = posY + zend;
			var posYzendzdem:int = posY + zend + zdem;
			
			var posZzdem:int = posZ + zdem;
			var posZzend:int = posZ + zend;
			var posZzendzdem:int = posZ + zend + zdem;
			
			var bin:FluidBIN;
			
			bin = binsYX;
			while(bin)
			{
				ind = bin.ind;
				Memory.writeFloat(
			    						Memory.readFloat(ind + posXzdem),
			    						ind + posX);
			    Memory.writeFloat(
			    						Memory.readFloat(ind + posXzend),
			    						ind + posXzendzdem);
			    						
			    Memory.writeFloat(
			    						Memory.readFloat(ind + posYzdem),
			    						ind + posY);
			    Memory.writeFloat(
			    						Memory.readFloat(ind + posYzend),
			    						ind + posYzendzdem);
			    						
				Memory.writeFloat(
			    						-Memory.readFloat(ind + posZzdem),
			    						ind + posZ);
			    Memory.writeFloat(
			    						-Memory.readFloat(ind + posZzend),
			    						ind + posZzendzdem);
				bin = bin.nextYX;        
			}
			
			var posXxdem:int = posX + xdem;
			var posXyend:int = posX + yend;
			var posXyendxdem:int = posX + yend + xdem;
			
			var posYxdem:int = posY + xdem;
			var posYyend:int = posY + yend;
			var posYyendxdem:int = posY + yend + xdem;
			
			var posZxdem:int = posZ + xdem;
			var posZyend:int = posZ + yend;
			var posZyendxdem:int = posZ + yend + xdem;
			
			bin = binsZX;
			while(bin)
			{
				ind = bin.ind;
				Memory.writeFloat(
			    						Memory.readFloat(ind + posXxdem),
			    						ind + posX);
			    Memory.writeFloat(
			    						Memory.readFloat(ind + posXyend),
			    						ind + posXyendxdem);
			    						
				Memory.writeFloat(
			    						-Memory.readFloat(ind + posYxdem),
			    						ind + posY);
			    Memory.writeFloat(
			    						-Memory.readFloat(ind + posYyend),
			    						ind + posYyendxdem);
			    						
				Memory.writeFloat(
			    						Memory.readFloat(ind + posZxdem),
			    						ind + posZ);
			    Memory.writeFloat(
			    						Memory.readFloat(ind + posZyend),
			    						ind + posZyendxdem);
				bin = bin.nextZX;        
			}
			
			var posX4:int = posX + 4;
			var posXxend:int = posX + xend;
			var posXxend4:int = posX + xend + 4;
			
			var posY4:int = posY + 4;
			var posYxend:int = posY + xend;
			var posYxend4:int = posY + xend + 4;
			
			var posZ4:int = posZ + 4;
			var posZxend:int = posZ + xend;
			var posZxend4:int = posZ + xend + 4;
			
			bin = binsZY;
			while(bin)
			{
				ind = bin.ind;
				Memory.writeFloat(
			    						-Memory.readFloat(ind + posX4),
			    						ind + posX);
			    Memory.writeFloat(
			    						-Memory.readFloat(ind + posXxend),
			    						ind + posXxend4);
			    
			    Memory.writeFloat(
			    						Memory.readFloat(ind + posY4),
			    						ind + posY);
			    Memory.writeFloat(
			    						Memory.readFloat(ind + posYxend),
			    						ind + posYxend4);
			    
			    Memory.writeFloat(
			    						Memory.readFloat(ind + posZ4),
			    						ind + posZ);
			    Memory.writeFloat(
			    						Memory.readFloat(ind + posZxend),
			    						ind + posZxend4);
				bin = bin.nextZY;        
			}
			
			Memory.writeFloat( 
								0.33 * (Memory.readFloat(posX4)
                                  + Memory.readFloat(posXxdem)
                                  + Memory.readFloat(posXzdem)),
                                  posX);
			Memory.writeFloat( 
								0.33 * (Memory.readFloat(posY4)
                                  + Memory.readFloat(posYxdem)
                                  + Memory.readFloat(posYzdem)),
                                  posY);
			Memory.writeFloat( 
								0.33 * (Memory.readFloat(posZ4)
                                  + Memory.readFloat(posZxdem)
                                  + Memory.readFloat(posZzdem)),
                                  posZ);
                                  
			Memory.writeFloat( 
								0.33 * (Memory.readFloat(posX+corner21)
                                  + Memory.readFloat(posXyend)
                                  + Memory.readFloat(posX + corner23)),
                                  posX + corner24);
			Memory.writeFloat( 
								0.33 * (Memory.readFloat(posY+corner21)
                                  + Memory.readFloat(posYyend)
                                  + Memory.readFloat(posY + corner23)),
                                  posY + corner24);
			Memory.writeFloat( 
								0.33 * (Memory.readFloat(posZ+corner21)
                                  + Memory.readFloat(posZyend)
                                  + Memory.readFloat(posZ + corner23)),
                                  posZ + corner24);
                                  
			Memory.writeFloat( 
								0.33 * (Memory.readFloat(posX+corner31)
                                  + Memory.readFloat(posX + corner32)
                                  + Memory.readFloat(posXzend)),//
                                  posX + corner34);
			Memory.writeFloat( 
								0.33 * (Memory.readFloat(posY+corner31)
                                  + Memory.readFloat(posY + corner32)
                                  + Memory.readFloat(posYzend)),//
                                  posY + corner34);
			Memory.writeFloat( 
								0.33 * (Memory.readFloat(posZ+corner31)
                                  + Memory.readFloat(posZ + corner32)
                                  + Memory.readFloat(posZzend)),//
                                  posZ + corner34);
                                  
			Memory.writeFloat( 
								0.33 * (Memory.readFloat(posX+corner41)//
                                  + Memory.readFloat(posX + corner42)
                                  + Memory.readFloat(posX + corner43)),
                                  posX + corner44);
			Memory.writeFloat( 
								0.33 * (Memory.readFloat(posY+corner41)//
                                  + Memory.readFloat(posY + corner42)
                                  + Memory.readFloat(posY + corner43)),
                                  posY + corner44);
			Memory.writeFloat( 
								0.33 * (Memory.readFloat(posZ+corner41)//
                                  + Memory.readFloat(posZ + corner42)
                                  + Memory.readFloat(posZ + corner43)),
                                  posZ + corner44);
                                  
			Memory.writeFloat( 
								0.33 * (Memory.readFloat(posXxend)
                                  + Memory.readFloat(posX + corner52)
                                  + Memory.readFloat(posX + corner53)),
                                  posX + corner54);
			Memory.writeFloat( 
								0.33 * (Memory.readFloat(posYxend)
                                  + Memory.readFloat(posY + corner52)
                                  + Memory.readFloat(posY + corner53)),
                                  posY + corner54);
			Memory.writeFloat( 
								0.33 * (Memory.readFloat(posZxend)
                                  + Memory.readFloat(posZ + corner52)
                                  + Memory.readFloat(posZ + corner53)),
                                  posZ + corner54);
                                  
			Memory.writeFloat( 
								0.33 * (Memory.readFloat(posX+corner61)
                                  + Memory.readFloat(posX + corner62)
                                  + Memory.readFloat(posX + corner63)),
                                  posX + corner64);
			Memory.writeFloat( 
								0.33 * (Memory.readFloat(posY+corner61)
                                  + Memory.readFloat(posY + corner62)
                                  + Memory.readFloat(posY + corner63)),
                                  posY + corner64);
			Memory.writeFloat( 
								0.33 * (Memory.readFloat(posZ+corner61)
                                  + Memory.readFloat(posZ + corner62)
                                  + Memory.readFloat(posZ + corner63)),
                                  posZ + corner64);
                                  
			Memory.writeFloat( 
								0.33 * (Memory.readFloat(posX+corner71)
                                  + Memory.readFloat(posX + corner72)
                                  + Memory.readFloat(posX + corner73)),
                                  posX + corner74);
			Memory.writeFloat( 
								0.33 * (Memory.readFloat(posY+corner71)
                                  + Memory.readFloat(posY + corner72)
                                  + Memory.readFloat(posY + corner73)),
                                  posY + corner74);
			Memory.writeFloat( 
								0.33 * (Memory.readFloat(posZ+corner71)
                                  + Memory.readFloat(posZ + corner72)
                                  + Memory.readFloat(posZ + corner73)),
                                  posZ + corner74);
                                  
			Memory.writeFloat( 
								0.33 * (Memory.readFloat(posX+corner81)
                                  + Memory.readFloat(posX + corner82)
                                  + Memory.readFloat(posX + corner83)),
                                  posX + corner84);
			Memory.writeFloat( 
								0.33 * (Memory.readFloat(posY+corner81)
                                  + Memory.readFloat(posY + corner82)
                                  + Memory.readFloat(posY + corner83)),
                                  posY + corner84);
			Memory.writeFloat( 
								0.33 * (Memory.readFloat(posZ+corner81)
                                  + Memory.readFloat(posZ + corner82)
                                  + Memory.readFloat(posZ + corner83)),
                                  posZ + corner84);
		}
		
		protected function buildBins():void
		{			
			var i:int, j:int, k:int;
			var ind:int;
			
			var p:FluidBIN;
			
			binIndex = new Vector.<FluidBIN>(bytesPerArray>>2, true);
		
			for (k = sizeZ; k > 0; --k) 
			{
		        for (j = sizeY; j > 0; --j) 
		        {
		            for (i = sizeX; i > 0; --i)
		            {
		            	ind = ((k * sizeXY2) + (j * sizeX2) + i) << 2;
		            	
		            	if( null == p )
		            	{
							p = bins = new FluidBIN();
		            	} else {
							p = p.next = new FluidBIN();
						}
						
						p.i = i;
						p.j = j;
						p.k = k;
						
						p.ind = ind;
						p.xi = posVelX + ind;
						p.yi = posVelY + ind;
						p.zi = posVelZ + ind;
						
						p.pxi = posPrevVelX + ind;
						p.pyi = posPrevVelY + ind;
						p.pzi = posPrevVelZ + ind;
						
						p.ri = posR + ind;
						p.gi = posG + ind;
						p.bi = posB + ind;
						
						p.pri = posPrevR + ind;
						p.pgi = posPrevG + ind;
						p.pbi = posPrevB + ind;
						
						binIndex[ind>>2] = p;
		            }
		        }
			}
			
			var yx:FluidBIN;
			
			for(j = sizeY; j > 0; --j) 
			{
			    for(i = sizeX; i > 0; --i) 
			    {
			    	ind = ((j * sizeX2) + i) << 2;

		    		if( null == yx )
	            	{
						yx = binsYX = new FluidBIN();
	            	} else {
						yx = yx.nextYX = new FluidBIN();
					}
					
					yx.ind = ind;
					yx.xi = posVelX + ind;
					yx.yi = posVelY + ind;
					yx.zi = posVelZ + ind;
					
					yx.pxi = posPrevVelX + ind;
					yx.pyi = posPrevVelY + ind;
					yx.pzi = posPrevVelZ + ind;
					
					yx.ri = posR + ind;
					yx.gi = posG + ind;
					yx.bi = posB + ind;
					
					yx.pri = posPrevR + ind;
					yx.pgi = posPrevG + ind;
					yx.pbi = posPrevB + ind;
					
					//binIndex[ind>>2] = yx;
				}
			}
			
			var zx:FluidBIN;
			
			for(k = sizeZ; k > 0; --k) 
			{
			    for(i = sizeX; i > 0; --i) 
			    {
			    	ind = ((k * sizeXY2) + i) << 2;

		    		if( null == zx )
	            	{
						zx = binsZX = new FluidBIN();
	            	} else {
						zx = zx.nextZX = new FluidBIN();
					}
					
					zx.ind = ind;
					zx.xi = posVelX + ind;
					zx.yi = posVelY + ind;
					zx.zi = posVelZ + ind;
					
					zx.pxi = posPrevVelX + ind;
					zx.pyi = posPrevVelY + ind;
					zx.pzi = posPrevVelZ + ind;
					
					zx.ri = posR + ind;
					zx.gi = posG + ind;
					zx.bi = posB + ind;
					
					zx.pri = posPrevR + ind;
					zx.pgi = posPrevG + ind;
					zx.pbi = posPrevB + ind;
					
					//binIndex[ind>>2] = zx;
				}
			}
			
			var zy:FluidBIN;
			
			for(k = sizeZ; k > 0; --k) 
			{
			    for(j = sizeY; j > 0; --j) 
			    {
			    	ind = ((k * sizeXY2) + j * sizeX2) << 2;

		    		if( null == zy )
	            	{
						zy = binsZY = new FluidBIN();
	            	} else {
						zy = zy.nextZY = new FluidBIN();
					}
					
					zy.ind = ind;
					zy.xi = posVelX + ind;
					zy.yi = posVelY + ind;
					zy.zi = posVelZ + ind;
					
					zy.pxi = posPrevVelX + ind;
					zy.pyi = posPrevVelY + ind;
					zy.pzi = posPrevVelZ + ind;
					
					zy.ri = posR + ind;
					zy.gi = posG + ind;
					zy.bi = posB + ind;
					
					zy.pri = posPrevR + ind;
					zy.pgi = posPrevG + ind;
					zy.pbi = posPrevB + ind;
					
					//binIndex[ind>>2] = zy;
				}
			}
			
			var fp:FluidBIN;
			
			for (k = sizeZ+1; k > -1; --k) 
			{
		        for (j = sizeY+1; j > -1; --j) 
		        {
		            for (i = sizeX+1; i > -1; --i)
		            {
		            	ind = ((k * sizeXY2) + (j * sizeX2) + i) << 2;
		            	
		            	if(binIndex[ind>>2]) continue;
		            	 
		            	if( null == fp )
		            	{
							fp = binsFull = new FluidBIN();
		            	} else {
							fp = fp.next = new FluidBIN();
						}
					
						fp.i = i;
						fp.j = j;
						fp.k = k;
						
						fp.ind = ind;
						fp.xi = posVelX + ind;
						fp.yi = posVelY + ind;
						fp.zi = posVelZ + ind;
						
						fp.pxi = posPrevVelX + ind;
						fp.pyi = posPrevVelY + ind;
						fp.pzi = posPrevVelZ + ind;
						
						fp.ri = posR + ind;
						fp.gi = posG + ind;
						fp.bi = posB + ind;
						
						fp.pri = posPrevR + ind;
						fp.pgi = posPrevG + ind;
						fp.pbi = posPrevB + ind;
						
						//binIndex[ind>>2] = fp;
		            }
		        }
			}
			fp = fp.next = bins;
			
			//
			initCorners();
		}
		
		protected function initCorners():void
		{
			var xdem:int = sizeX2 << 2;
			var xend:int = sizeX << 2;
			var yend:int = xdem * sizeY;
			var zdem:int = sizeXY2 << 2;
			var zend:int = zdem * sizeZ;
			
			corner21 =  4+yend+xdem;
            corner23 =  yend+xdem + zdem;
            corner24 =  yend+xdem;
            
			corner31 =  4+zend+zdem;
            corner32 =  xdem + zend+zdem;
            corner34 =  zend + zdem;
            
			corner41 =  4 + zend + zdem + yend+xdem;
            corner42 =  yend + zend + zdem;
            corner43 =  yend+xdem+zend;
            corner44 =  zend + zdem + yend+xdem;
            
            corner52 =  xend+4 + xdem;
            corner53 =  xend+4 + zdem;
            corner54 =  xend+4;
            
			corner61 =  xend + yend+xdem;
            corner62 =  xend+4 + yend;
            corner63 =  xend+4 + yend+xdem + zdem;
            corner64 =  xend+4 + yend+xdem;
            
			corner71 =  xend + zend+zdem;
            corner72 =  xend+4 + xdem + zend+zdem;
            corner73 =  xend+4 + zend;
            corner74 =  xend+4 + zend+zdem;
            
			corner81 =  xend + yend+xdem + zend+zdem;
            corner82 =  xend+4 + yend + zend+zdem;
            corner83 =  xend+4 + yend+xdem + zend;
            corner84 =  xend+4 + zend+zdem + yend+xdem;
		}
		
		protected var corner21:int;
		protected var corner23:int;
		protected var corner24:int;
		
		protected var corner31:int;
		protected var corner32:int;
		protected var corner34:int;
		
		protected var corner41:int;
		protected var corner42:int;
		protected var corner43:int;
		protected var corner44:int;
		
		protected var corner52:int;
		protected var corner53:int;
		protected var corner54:int;
		
		protected var corner61:int;
		protected var corner62:int;
		protected var corner63:int;
		protected var corner64:int;
		
		protected var corner71:int;
		protected var corner72:int;
		protected var corner73:int;
		protected var corner74:int;
		
		protected var corner81:int;
		protected var corner82:int;
		protected var corner83:int;
		protected var corner84:int;
	}
}

internal final class FluidBIN
{
	public var next:FluidBIN;
	public var nextYX:FluidBIN;
	public var nextZX:FluidBIN;
	public var nextZY:FluidBIN;
	
	public var i:int;
	public var j:int;
	public var k:int;
	
	public var ind:int;
	
	public var xi:int;
	public var yi:int;
	public var zi:int;
	
	public var pxi:int;
	public var pyi:int;
	public var pzi:int;
	
	public var ri:int;
	public var gi:int;
	public var bi:int;
	
	public var pri:int;
	public var pgi:int;
	public var pbi:int;
}