package ru.inspirit.steering.behavior 
{
	import ru.inspirit.steering.Vehicle;
	import ru.inspirit.steering.behavior.AbstractBehavior;

	/**
	 * @author Eugene Zatepyakin
	 */
	public final class UnalignedCollisionAvoidance extends AbstractBehavior 
	{
		public var collisionDangerThreshold:Number = 9; // should be around 3 * Vehicle Radius
		public var predict:Number = 25;
		
		public var avoidList:Vehicle;
		
		public function UnalignedCollisionAvoidance(avoidList:Vehicle = null, collisionDangerThreshold:Number = 9, predict:Number = 25)
		{
			this.avoidList = avoidList;
			this.collisionDangerThreshold = collisionDangerThreshold;
			this.predict = predict;
		}
		
		override public function apply(veh:Vehicle):void
		{			
			var a:Number, b:Number, c:Number, t:Number;
			var dx:Number, dy:Number, dz:Number;
			var pdx:Number, pdy:Number, pdz:Number;
			var vx:Number = veh.velocity.x;
			var vy:Number = veh.velocity.y;
			var vz:Number = veh.velocity.z;
			var px:Number = veh.position.x;
			var py:Number = veh.position.y;
			var pz:Number = veh.position.z;
			var mag:Number;
			var d:Number;
			
			var threat:Vehicle = null;
			var time:Number;
			var minTime:Number = 20;
			var relSpeed:Number;
			var other:Vehicle = avoidList;
			
			while(other) 
			{
				if(other == veh)
				{
					other = other.next;
					continue;
				}
				// predict Nearest Approach Time
				dx = other.velocity.x - vx;
				dy = other.velocity.y - vy;
				dz = other.velocity.z - vz;
				
				a = dx;
				if (a < 0) a = -a;
				
				b = dy;
				if (b < 0) b = -b;
				
				c = dz;
				if (c < 0) c = -c;
				
				if(a < b) {
					t = a;
					a = b;
					b = t;
				}
				if(a < c) {
					t = a;
					a = c;
					c = t;
				}
				
				relSpeed = a * 0.9375 + (b + c) * 0.375;
				
				mag = Math.sqrt( dx * dx + dy * dy + dz * dz );
				if(mag != 0) 
				{
					mag = 1 / mag;
					dx *= mag;
					dy *= mag;
					dz *= mag;
				}
				
				pdx = px - other.position.x;
				pdy = py - other.position.y;
				pdz = pz - other.position.z;
				
				time = (dx * pdx + dy * pdy + dz * pdz) / relSpeed;
				
				if(time >= 0 && time < minTime)
				{
					// compute Nearest Approach Positions
					dx = (vx * time + px) - (other.velocity.x * time + other.position.x);
					dy = (vy * time + py) - (other.velocity.y * time + other.position.y);
					dz = (vz * time + pz) - (other.velocity.z * time + other.position.z);
					
					a = dx;
					if (a < 0) a = -a;
					
					b = dy;
					if (b < 0) b = -b;
					
					c = dz;
					if (c < 0) c = -c;
					
					if(a < b) {
						t = a;
						a = b;
						b = t;
					}
					if(a < c) {
						t = a;
						a = c;
						c = t;
					}
					
					d = a * 0.9375 + (b + c) * 0.375;
					
					if(d < collisionDangerThreshold)
					{
						minTime = time;
						threat = other;
					}
				}
				
				other = other.next;
			}
			
			if(threat)
			{
				pdx = (threat.velocity.x * minTime + threat.position.x);
				pdy = (threat.velocity.y * minTime + threat.position.y);
				pdz = (threat.velocity.z * minTime + threat.position.z);
				dx = px - pdx;
				dy = py - pdy;
				dz = pz - pdz;
				d = (dx * veh.up.x + dy * veh.up.y + dz * veh.up.z);
				dx = d * veh.up.x;
				dy = d * veh.up.y;
				dz = d * veh.up.z;
				
				if(veh.velocity.approximateLength() < 0.2 * veh.maxSpeed)
				{
					accumulator.x = dx + veh.forward.x;
					accumulator.y = dy + veh.forward.y;
					accumulator.z = dz + veh.forward.z;
				} else {
					pdx = (vx * minTime + px) - pdx;
					pdy = (vy * minTime + py) - pdy;
					pdz = (vz * minTime + pz) - pdz;
					d = (pdx * veh.forward.x + pdy * veh.forward.y + pdz * veh.forward.z);
					
					accumulator.x = dx + veh.forward.x * d;
					accumulator.y = dy + veh.forward.y * d;
					accumulator.z = dz + veh.forward.z * d;
				}
				veh.applyGlobalForce(accumulator);
			}
		}
	}
}
