package ru.inspirit.fluid.emitter 
{
	import ru.inspirit.fluid.FluidSolver3D;
	import ru.inspirit.fluid.Particle3DManager;

	import com.joa_ebert.apparat.memory.Memory;

	/**
	 * @author Eugene Zatepyakin
	 */
	public class FluidForceEmitter implements IFluidEmitter
	{
		protected const RAD:Number = Math.PI / 180;
		
		public var x:int;
		public var y:int;
		public var z:int;
		
		public var fx:Number;
		public var fy:Number;
		public var fz:Number;
		
		public var solver:FluidSolver3D;
		public var particles:Particle3DManager;
		
		public var index:int;
		
		public var rotXLight:Number = Math.random()*360;
		public var rotYLight:Number = Math.random()*360;
		
		protected var posVelX:int;
		protected var posVelY:int;
		protected var posVelZ:int;
		
		public var forceMultiplier:Number = 0.5;
		
		public function FluidForceEmitter(solver:FluidSolver3D, particles:Particle3DManager = null, x:int = 0, y:int = 0, z:int = 0)
		{
			this.solver = solver;
			this.particles = particles;
			
			this.x = x;
			this.y = y;
			this.z = z;
			
			this.index = ((z * (solver.sizeX + 2) * (solver.sizeY + 2)) + (y * (solver.sizeX + 2)) + x) << 2;
			
			posVelX = solver.posVelX + index;
			posVelY = solver.posVelY + index;
			posVelZ = solver.posVelZ + index;
		}
		
		public function update():void
		{			
			rotXLight += 1;
			rotYLight += 3;
			var radX:Number = rotXLight * RAD;
			var radY:Number = rotYLight * RAD;
			
			fx = Math.cos(radX) * Math.cos(radY) * forceMultiplier;
			fy = -Math.sin(radY) * forceMultiplier;
			fz = Math.sin(radX) * Math.cos(radY) * forceMultiplier;
			
			particles.addParticles(x, y, z, 10);
			
			Memory.writeFloat(
								Memory.readFloat( posVelX ) + fx, 
								posVelX);
			Memory.writeFloat(
								Memory.readFloat( posVelY ) + fy, 
								posVelY);
			Memory.writeFloat(
								Memory.readFloat( posVelZ ) + fz, 
								posVelZ);
		}
		
		public function updateWithForce(fx:Number, fy:Number, fz:Number):void
		{
			//
		}
		public function setRandomPosition():void
		{
			
		}
	}
}
