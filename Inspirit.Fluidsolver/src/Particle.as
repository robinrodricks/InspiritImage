package  
{
	import ru.inspirit.utils.Random;
	import ru.inspirit.utils.FluidSolver;
	
	/**
	 * Basic particle unit
	 * @author Eugene Zatepyakin
	 */
	public class Particle 
	{
		public static var MOMENTUM:Number = 0.5;
		public static var FLUID_FORCE:Number = 0.6;
		
		private var fs:FluidSolver;
		private var sw:uint;
		private var sh:uint;
		private var isw:Number;
		private var ish:Number;
		
		public var x:Number, y:Number;
		public var vx:Number, vy:Number;
		public var radius:Number;
		public var alpha:Number;
		public var mass:Number;
		
		public var next:Particle;
		
		public function Particle() 
		{
			fs = Main.fSolver;
			sw = Main.sw;
			sh = Main.sh;
			isw = Main.isw;
			ish = Main.ish;
			alpha = 0;
		}
		
		public function init(x:Number = 0, y:Number = 0):void
		{
			this.x = x;
			this.y = y;
			vx = vy = 0;
			radius = 5;
			alpha = Random.float(.3, 1);
			mass = Random.float(.1, 1);
		}
		
		public function update():void
		{
			if (alpha == 0) return;
			
			const fluidIndex:int = fs.getIndexForNormalizedPosition(x * isw, y * ish);
			
			vx = fs.u[fluidIndex] * sw * mass * FLUID_FORCE + vx * MOMENTUM;
			vy = fs.v[fluidIndex] * sh * mass * FLUID_FORCE + vy * MOMENTUM;
			
			x += vx;
			y += vy;
			
			if (x < 0) {
				if (fs.wrapX) {
					x += sw;
				} else {
					x = 1;
					vx *= -1;
				}
			}
			else if (x > sw) {
				if (fs.wrapX) {
					x -= sw;
				} else {
					x = sw - 1;
					vx *= -1;
				}
			}
			
			if (y < 0) {
				if (fs.wrapY) {
					y += sh;
				} else {
					y = 1;
					vy *= -1;
				}
			}
			else if (y > sh) {
				if (fs.wrapY) {
					y -= sh;
				} else {
					y = sh - 1;
					vy *= -1;
				}
			}
			
			// hackish way to make particles glitter when the slow down a lot
			if(vx * vx + vy * vy < .5) {
				alpha = 0;
				return;
				vx = Random.float(-1, 1);
				vy = Random.float(-1, 1);
			}
			
			// fade out a bit (and kill if alpha == 0);
			alpha *= 0.999;
			if(alpha < 0.01) alpha = 0;
		}
		
	}
	
}