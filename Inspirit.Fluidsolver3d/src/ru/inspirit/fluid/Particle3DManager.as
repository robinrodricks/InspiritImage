package ru.inspirit.fluid 
{
	import flash.geom.Matrix3D;
	import flash.display.BitmapData;
	import com.joa_ebert.apparat.memory.Memory;

	/**
	 * @author Eugene Zatepyakin
	 */
	public class Particle3DManager 
	{
		public static const MAX_PARTICLES:int = 10000;
		public static const MOMENTUM:Number = 0.5;
		public static const FLUID_FORCE:Number = 0.6;
		
		public static var seed:uint = 1;//int(Math.random() * 0x7FFFFFFE) + 1;
		
		public var head:Particle3D;
		public var tail:Particle3D;
		
		public var numUsedParticles:int = 0;
		
		public var solver:FluidSolver3D;
		
		public var xMin:Number;
		public var yMin:Number;
		public var zMin:Number;
		public var xMax:Number;
		public var yMax:Number;
		public var zMax:Number;
		
		protected var invW:Number;
		protected var invH:Number;
		protected var invD:Number;
		
		protected var invWfsW:Number;
		protected var invHfsH:Number;
		protected var invDfsD:Number;
		
		protected var sW:int;
		protected var sH:int;
		protected var sD:int;
		protected var hsW:int;
		protected var hsH:int;
		protected var hsD:int;
		
		protected var fsW:int;
		protected var fsH:int;
		protected var fsD:int;
		protected var fsW2:int;
		protected var fsWfsH2:int;
		
		protected var posVelX:int;
		protected var posVelY:int;
		protected var posVelZ:int;
		
		public function Particle3DManager(solver:FluidSolver3D, xMin:Number = 0, yMin:Number = 0, xMax:Number = 800, yMax:Number = 600, zMin:Number = 0, zMax:Number = 300)
		{
			head = tail = new Particle3D();
			
			this.solver = solver;
			
			fsW = solver.sizeX;
			fsH = solver.sizeY;
			fsD = solver.sizeZ;
			
			fsW2 = fsW + 2;
			fsWfsH2 = fsW2 * (fsH+2);
			
			sW = xMax - xMin;
			sH = yMax - yMin;
			sD = zMax - zMin;
			
			invW = 1 / sW;
			invH = 1 / sH;
			invD = 1 / sD;
			
			hsW = sW * 0.5;
			hsH = sH * 0.5;
			hsD = sD * 0.5;
			
			this.xMin = -hsW;
			this.xMax = hsW;
			this.yMin = -hsH;
			this.yMax = hsH;
			this.zMin = -hsD;
			this.zMax = hsD;
			
			posVelX = solver.posVelX;
			posVelY = solver.posVelY;
			posVelZ = solver.posVelZ;
			
			invWfsW = invW*fsW;
			invHfsH = invH*fsH;
			invDfsD = invD*fsD;
		}

		public final function render(bmp:BitmapData, matrix:Matrix3D, focalLength:Number = 60):void
		{
			var p:Particle3D = head.next;
			var pp:Particle3D = head;
			var vx:Number, vy:Number, vz:Number;
			var x:Number, y:Number, z:Number, w:Number;
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
			
			var sw:int = bmp.width;
			var sh:int = bmp.height;
			var cx: Number = sw * 0.5;
			var cy: Number = sh * 0.5;
			
			var i:int, j:int, k:int, a:int;
			var c:uint;
			var ind:int;
			
			
			while(p)
			{
				if(p.alpha > 0.1)
				{
					i = 1 + (p.x + hsW)*invWfsW;
					j = 1 + (p.y + hsH)*invHfsH;
					k = 1 + (p.z + hsD)*invDfsD;
					
					ind = ((k * fsWfsH2) + (j * fsW2) + i) << 2;
					
					vx = Memory.readFloat(posVelX + ind) * sW * p.mass * FLUID_FORCE + p.vx * MOMENTUM;
					vy = Memory.readFloat(posVelY + ind) * sH * p.mass * FLUID_FORCE + p.vy * MOMENTUM;
					vz = Memory.readFloat(posVelZ + ind) * sD * p.mass * FLUID_FORCE + p.vz * MOMENTUM;
					
					x = p.x + vx;
					y = p.y + vy;
					z = p.z + vz;
					
					if(x < xMin)
					{
						x = xMin + 2;
						vx *= -1;
					} else if(x > xMax)
					{
						x = xMax - 2;
						vx *= -1;
					}
					
					if(y < yMin)
					{
						y = yMin + 2;
						vy *= -0.5;
					} else if(y > yMax)
					{
						y = yMax - 2;
						vy *= -1;
					}
					
					if(z < zMin)
					{
						z = zMin + 2;
						vz *= -1;
					} else if(z > zMax)
					{
						z = zMax - 2;
						vz *= -1;
					}
					
					pz = focalLength + x * p02 + y * p12 + z * p22 + p32;
					
					//if(pz > 0)
					//{
						i = int( ( w = focalLength / pz ) * ( x * p00 + y * p10 + z * p20) + cx );
						j = int( w * ( x * p01 + y * p11 + z * p21) + cy );
						
						a = p.alpha * 0xFF;
						c = a << 24 | a << 16 | a << 8 | a;
						
						if(i > -1 && i < sw && j > -1 && j < sh) 
						{
							// inlined line drawing
							var vxx:Number = x - vx;
							var vyy:Number = y - vy;
							var vzz:Number = z - vz;
							if(vxx < xMin) vxx = xMin;
							if(vxx > xMax) vxx = xMax;
							if(vyy < yMin) vyy = yMin;
							if(vyy > yMax) vyy = yMax;
							if(vzz < zMin) vzz = zMin;
							if(vzz > zMax) vzz = zMax;
							var x1:int = int( w * ( (vxx) * p00 + (vyy) * p10 + (vzz) * p20) + cx );
							var y1:int = int( w * ( (vxx) * p01 + (vyy) * p11 + (vzz) * p21) + cy );
							
							var shortLen:int = j - y1;
							var longLen:int = i - x1;
							var inc:int;
							var multDiff:Number;
						
							if((shortLen ^ (shortLen >> 31)) - (shortLen >> 31) > (longLen ^ (longLen >> 31)) - (longLen >> 31))
							{
								shortLen ^= longLen;
								longLen ^= shortLen;
								shortLen ^= longLen;
						
								inc = longLen < 0 ? -1 : 1;
								multDiff = longLen == 0 ? shortLen : shortLen / longLen;
								
								for (i = 0; i != longLen; i += inc)
								{
									bmp.setPixel32(x1 + i*multDiff, y1+i, c);
							    }
							}
							else
							{
								inc = longLen < 0 ? -1 : 1;
								multDiff = longLen == 0 ? shortLen : shortLen / longLen;
								
								for (i = 0; i != longLen; i += inc)
								{
									bmp.setPixel32(x1+i, y1+i*multDiff, c);
								}
							}
						}
					//}
					
					p.x = x;
					p.y = y;
					p.z = z;
					
					p.vx = vx;
					p.vy = vy;
					p.vz = vz;
					
					p.alpha *= 0.996;
					
					pp = p;
					p = p.next;
				} 
				else 
				{
					pp.next = p.next;
					p.next = null;
					p = null;
					p = pp.next;
					numUsedParticles--;
				}
			}
		}
		
		public final function addParticles(x:Number, y:Number, z:Number, n:int):void
		{
			var p:Particle3D;
			
			while(--n > -1)
			{
				if(numUsedParticles == MAX_PARTICLES) break;
				
				p = new Particle3D();
				numUsedParticles++;
				
				p.x = x + random()*20 - 10;
				p.y = y + random()*20 - 10;
				p.z = z + random()*20 - 10;
				
				p.alpha = random()*0.7 + 0.3;
				p.mass = random()*0.9 + 0.1;
				
				tail = tail.next = p;
			}
		}
		
		public final function removeParticle(p:Particle3D):void
		{
			var node:Particle3D = head;
			while (node.next != p) node = node.next;
            //if (node.next == tail) tail = node;
            node.next = p.next;
		}

		public final function line(bmp:BitmapData, x1:int, y1:int, x2:int, y2:int, c:uint):void
		{	
			var i:int;
			var shortLen:int = y2 - y1;
			var longLen:int = x2 - x1;
		
			if((shortLen ^ (shortLen >> 31)) - (shortLen >> 31) > (longLen ^ (longLen >> 31)) - (longLen >> 31))
			{
				shortLen ^= longLen;
				longLen ^= shortLen;
				shortLen ^= longLen;
		
				var yLonger:Boolean = true;
			}
			else
			{
				yLonger = false;
			}
		
			var inc:int = longLen < 0 ? -1 : 1;
			var multDiff:Number = longLen == 0 ? shortLen : shortLen / longLen;
		
			if (yLonger)
			{
				for (i = 0; i != longLen; i += inc)
				{
					bmp.setPixel32(x1 + i*multDiff, y1+i, c);
			    }
			}
			else
			{
				for (i = 0; i != longLen; i += inc)
				{
					bmp.setPixel32(x1+i, y1+i*multDiff, c);
				}
			}
		}
		
		public static function random():Number
        {
        	return (seed = (seed * 16807) % 2147483647) / 2147483647;
        }
	}
}