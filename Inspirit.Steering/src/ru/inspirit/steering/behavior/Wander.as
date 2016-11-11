package ru.inspirit.steering.behavior 
{
	import ru.inspirit.steering.SteerVector3D;
	import ru.inspirit.steering.Vehicle;
	import ru.inspirit.steering.behavior.AbstractBehavior;

	/**
	 * @author Eugene Zatepyakin
	 */
	public final class Wander extends AbstractBehavior 
	{
		protected const SQRT2:Number = Math.SQRT2;
		
		protected var wanderJiggle:SteerVector3D = new SteerVector3D();
		protected var wanderGlobal:SteerVector3D = new SteerVector3D();
		
		public var wanderStrength:Number = .2;
		public var wanderRate:Number = 0.8;
		public var wanderDirection:SteerVector3D = new SteerVector3D();
		
		public function Wander(wanderStrength:Number = 0.2, wanderRate:Number = 0.8)
		{
			this.wanderStrength = wanderStrength;
			this.wanderRate = wanderRate;
		}
		
		override public function apply(veh:Vehicle):void
		{
			wanderJiggle.setUnitRandom();
			
			wanderJiggle.x *= wanderRate;
			wanderJiggle.y *= wanderRate;
			wanderJiggle.z *= wanderRate;
			
			//wanderDirection.x += wanderJiggle.x;
			wanderDirection.y += wanderJiggle.y;
			wanderDirection.z += wanderJiggle.z;
			
			wanderDirection.setApproximateNormalize();
			
			veh.globalizeDirection(wanderDirection, wanderGlobal);
			
			wanderGlobal.x *= wanderStrength;
			wanderGlobal.y *= wanderStrength;
			wanderGlobal.z *= wanderStrength;
			
			accumulator.x = veh.forward.x * SQRT2 + wanderGlobal.x;
			accumulator.y = veh.forward.y * SQRT2 + wanderGlobal.y;
			accumulator.z = veh.forward.z * SQRT2 + wanderGlobal.z;
			
			veh.applyGlobalForce(accumulator);
		}
	}
}
