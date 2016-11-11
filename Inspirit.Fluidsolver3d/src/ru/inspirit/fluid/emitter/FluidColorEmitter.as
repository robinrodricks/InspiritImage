package ru.inspirit.fluid.emitter 
{
	import ru.inspirit.fluid.FluidSolver3D;
	import ru.inspirit.utils.ColorUtils;

	import com.joa_ebert.apparat.memory.Memory;

	/**
	 * @author Eugene Zatepyakin
	 */
	public class FluidColorEmitter implements IFluidEmitter
	{
		public var x:int;
		public var y:int;
		public var z:int;
		
		public var solver:FluidSolver3D;
		
		public var index:int;
		public var frameCount:int = 0;
		
		protected var nx:Number;
		protected var ny:Number;
		protected var nz:Number;
		
		protected var posR:int;
		protected var posG:int;
		protected var posB:int;
		
		public function FluidColorEmitter(solver:FluidSolver3D, x:int = 0, y:int = 0, z:int = 0)
		{
			this.solver = solver;
			
			this.x = x;
			this.y = y;
			this.z = z;
			
			this.nx = x / solver.sizeX;
			this.ny = y / solver.sizeY;
			this.nz = z / solver.sizeZ;
			
			this.index = ((z * (solver.sizeX + 2) * (solver.sizeY + 2)) + (y * (solver.sizeX + 2)) + x) << 2;
			
			posR = solver.posR + index;
			posG = solver.posG + index;
			posB = solver.posB + index;
		}
		
		public function update():void
		{			
			const colorMult:Number = 30;
			const hue:Number = ((nx + ny + nz) * 120 + frameCount) % 360;
			const rgb:Object = ColorUtils.HSB2GRB(hue);
			
			Memory.writeFloat(
								Memory.readFloat( posR ) + rgb.r * colorMult, 
								posR);
			Memory.writeFloat(
								Memory.readFloat( posG ) + rgb.g * colorMult, 
								posG);
			Memory.writeFloat(
								Memory.readFloat( posB ) + rgb.b * colorMult, 
								posB);
			
			frameCount++;
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
