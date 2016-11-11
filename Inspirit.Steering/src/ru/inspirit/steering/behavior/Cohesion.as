package ru.inspirit.steering.behavior 
{
	import ru.inspirit.steering.Vehicle;
	import ru.inspirit.steering.behavior.AbstractBehavior;

	/**
	 * @author Eugene Zatepyakin
	 */
	public final class Cohesion extends AbstractBehavior 
	{
		public var cohereDist:Number;
		public var cohereAngleCos:Number;
		public var cohereStrength:Number;
		
		public var cohereList:Vehicle;
		
		public function Cohesion(cohereList:Vehicle = null, cohereDist:Number = 15, cohereAngle:Number = 150, cohereStrength:Number = 1)
		{
			this.cohereList = cohereList;
			this.cohereDist = cohereDist;
			this.cohereAngleCos = Math.cos(((cohereAngle / 360) * 2.0) * Math.PI);
			this.cohereStrength = cohereStrength;
		}
		
		override public function apply(veh:Vehicle):void
		{
			var ax:Number = 0, ay:Number = 0, az:Number = 0;
			
			var count:int = 0;
			
			var other:Vehicle = cohereList;
			
			while(other) 
			{
				if(other == veh)
				{
					other = other.next;
					continue;
				}
				
				if( nearby(veh, other, cohereDist, cohereAngleCos) )
				{
					count++;
					
					ax += other.position.x;
					ay += other.position.y;
					az += other.position.z;
				}
				other = other.next;
			}
			
			if(count)
			{
				var f:Number = 1 / count;
				ax *= f;
				ay *= f;
				az *= f;
				
				accumulator.x = ax - veh.position.x;
				accumulator.y = ay - veh.position.y;
				accumulator.z = az - veh.position.z;
				accumulator.setApproximateNormalize();
			
				accumulator.x *= cohereStrength;
				accumulator.y *= cohereStrength;
				accumulator.z *= cohereStrength;
				
				veh.applyGlobalForce(accumulator);
			}
		}
	}
}
