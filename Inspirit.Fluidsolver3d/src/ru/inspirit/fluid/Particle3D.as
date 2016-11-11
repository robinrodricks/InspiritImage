package ru.inspirit.fluid 
{

	/**
	 * @author Eugene Zatepyakin
	 */
	public final class Particle3D 
	{		
		public var next:Particle3D = null;
		
		public var x:Number;
		public var y:Number;
		public var z:Number;
		
		public var vx:Number = 0;
		public var vy:Number = 0;
		public var vz:Number = 0;
		
		public var alpha:Number;
		public var mass:Number;
	}
}
