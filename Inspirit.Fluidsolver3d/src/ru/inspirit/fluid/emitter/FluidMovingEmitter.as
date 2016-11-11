package ru.inspirit.fluid.emitter 
{
	import ru.inspirit.fluid.FluidSolver3D;
	import ru.inspirit.fluid.Particle3DManager;
	import ru.inspirit.utils.ColorUtils;

	import com.joa_ebert.apparat.memory.Memory;

	/**
	 * @author Eugene Zatepyakin
	 */
	public class FluidMovingEmitter implements IFluidEmitter
	{
		protected const RAD:Number = Math.PI / 180;
		
		public var x:int;
		public var y:int;
		public var z:int;
		
		public var fx:Number;
		public var fy:Number;
		public var fz:Number;
		
		protected var nx:Number;
		protected var ny:Number;
		protected var nz:Number;
		
		protected var rx:Number;
		protected var ry:Number;
		protected var rz:Number;
		protected var dx:Number = 1;
		protected var dy:Number = -1;
		protected var dz:Number = 1;
		
		public var solver:FluidSolver3D;
		public var particles:Particle3DManager;
		
		public var index:int;
		
		public var rotXLight:Number = Math.random()*360;
		public var rotYLight:Number = Math.random()*360;
		
		public var frameCount:int = 0;
		
		protected var posVelX:int;
		protected var posVelY:int;
		protected var posVelZ:int;
		
		protected var posR:int;
		protected var posG:int;
		protected var posB:int;
		
		protected var sizeX:int;
		protected var sizeY:int;
		protected var sizeZ:int;
		
		protected var sizeX2:int;
		protected var sizeY2:int;
		protected var sizeXY2:int;
		
		protected var hsizeX:int;
		protected var hsizeY:int;
		protected var hsizeZ:int;
		
		public function FluidMovingEmitter(solver:FluidSolver3D, particles:Particle3DManager = null, x:int = 0, y:int = 0, z:int = 0)
		{
			this.solver = solver;
			this.particles = particles;
			
			posVelX = solver.posVelX;
			posVelY = solver.posVelY;
			posVelZ = solver.posVelZ;
			
			posR = solver.posR;
			posG = solver.posG;
			posB = solver.posB;
			
			sizeX = solver.sizeX;
			sizeY = solver.sizeY;
			sizeZ = solver.sizeZ;
			
			hsizeX = solver.sizeX >> 1;
			hsizeY = solver.sizeY >> 1;
			hsizeZ = solver.sizeZ >> 1;
			
			sizeX2 = solver.sizeX+2;
			sizeY2 = solver.sizeY+2;
			
			sizeXY2 = sizeX2 * sizeY2;
			
			this.x = x;
			this.y = y;
			this.z = z;
		}
		
		public function update():void
		{
			rotXLight += 1;
			rotYLight += 3;
			var radX:Number = rotXLight * RAD;
			var radY:Number = rotYLight * RAD;
			
			var px:Number = x;
			var py:Number = y;
			var pz:Number = z;
			
			x = (Math.cos(radX) * Math.cos(radY) * hsizeX) + hsizeX;
			y = (-Math.sin(radY) * hsizeY) + hsizeY;
			z = (Math.sin(radX) * Math.cos(radY) * hsizeZ) + hsizeZ;
			if(x < 1) x = 1;
			if(y < 1) y = 1;
			if(z < 1) z = 1;
			if(x > sizeX) x = sizeX;
			if(y > sizeY) y = sizeY;
			if(z > sizeZ) z = sizeZ;
			
			index = ((z * sizeXY2) + (y * sizeX2) + x) << 2;
			
			this.fx = (x - px) * 0.8;
			this.fy = (y - py) * 0.8;
			this.fz = (z - pz) * 0.8;
			
			var pos:int;
			Memory.writeFloat(
								Memory.readFloat( (pos = posVelX + index) ) + fx, 
								pos);
			Memory.writeFloat(
								Memory.readFloat( (pos = posVelY + index) ) + fy, 
								pos);
			Memory.writeFloat(
								Memory.readFloat( (pos = posVelZ + index) ) + fz, 
								pos);
								
			nx = x / sizeX;
			ny = y / sizeY;
			nz = z / sizeZ;
			
			if(frameCount & 1) particles.addParticles(nx*460-230, ny*460-230, nz*460-230, 15);
			
			const colorMult:Number = 30;
			const hue:Number = ((nx + ny + nz) * 120 + frameCount) % 360;
			const rgb:Object = ColorUtils.HSB2GRB(hue);
			
			Memory.writeFloat(
								Memory.readFloat( (pos = posR + index) ) + rgb.r * colorMult, 
								pos);
			Memory.writeFloat(
								Memory.readFloat( (pos = posG + index) ) + rgb.g * colorMult, 
								pos);
			Memory.writeFloat(
								Memory.readFloat( (pos = posB + index) ) + rgb.b * colorMult, 
								pos);
			
			frameCount++;
			
		}
		
		public function updateWithForce(tx:Number, ty:Number, tz:Number):void
		{
			var mult:Number = 1;
			
			var px:Number = x;
			var py:Number = y;
			var pz:Number = z;
			
			/*rx += dx * tx * mult;
			ry += dy * ty * mult;
			rz += dz * tz * mult;
			
			if(rx < 1) {rx = 1; dx *= -1;}
			if(ry < 1) {ry = 1; dy *= -1;}
			if(rz < 1) {rz = 1; dz *= -1;}
			if(rx > sizeX) {rx = sizeX; dx *= -1;}
			if(ry > sizeY) {ry = sizeY; dy *= -1;}
			if(rz > sizeZ) {rz = sizeZ; dz *= -1;}*/
			
			x = tx * sizeX;
			y = ty * sizeY;
			z = tz * sizeZ;
			//x = rx;
			//y = ry;
			//z = rz;
			if(x < 1) x = 1;
			if(y < 1) y = 1;
			if(z < 1) z = 1;
			if(x > sizeX) x = sizeX;
			if(y > sizeY) y = sizeY;
			if(z > sizeZ) z = sizeZ;
			
			index = ((z * sizeXY2) + (y * sizeX2) + x) << 2;
			mult = 0.8;
			this.fx = (x - px) * mult;
			this.fy = (y - py) * mult;
			this.fz = (z - pz) * mult;
			
			var pos:int;
			Memory.writeFloat(
								Memory.readFloat( (pos = posVelX + index) ) + fx, 
								pos);
			Memory.writeFloat(
								Memory.readFloat( (pos = posVelY + index) ) + fy, 
								pos);
			Memory.writeFloat(
								Memory.readFloat( (pos = posVelZ + index) ) + fz, 
								pos);
								
			nx = x / sizeX;
			ny = y / sizeY;
			nz = z / sizeZ;
			
			if(frameCount & 1) particles.addParticles(nx*460-230, ny*460-230, nz*460-230, 15);
			
			const colorMult:Number = 30;
			const hue:Number = ((nx + ny + nz) * 120 + frameCount) % 360;
			const rgb:Object = ColorUtils.HSB2GRB(hue);
			
			Memory.writeFloat(
								Memory.readFloat( (pos = posR + index) ) + rgb.r * colorMult, 
								pos);
			Memory.writeFloat(
								Memory.readFloat( (pos = posG + index) ) + rgb.g * colorMult, 
								pos);
			Memory.writeFloat(
								Memory.readFloat( (pos = posB + index) ) + rgb.b * colorMult, 
								pos);
			
			frameCount++;
		}
		
		public function setRandomPosition():void
		{
			x = 1 + Math.random() * sizeX;
			y = 1 + Math.random() * sizeY;
			z = 1 + Math.random() * sizeZ;
			
			rx = x;// sizeX;
			ry = y;// sizeY;
			rz = z;// sizeZ;
		}
	}
}
