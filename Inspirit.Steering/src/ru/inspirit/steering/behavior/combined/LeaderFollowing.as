package ru.inspirit.steering.behavior.combined 
{
	import ru.inspirit.steering.SteerVector3D;
	import ru.inspirit.steering.Vehicle;
	import ru.inspirit.steering.behavior.AbstractBehavior;

	/**
	 * Leader following bahevior
	 * Should update to allow users control leader vehicle
	 * 
	 * @author Eugene Zatepyakin
	 */
	public final class LeaderFollowing extends AbstractBehavior 
	{
		// WONDER OPTIONS
		protected const SQRT2:Number = Math.SQRT2;
		
		protected var wanderJiggle:SteerVector3D = new SteerVector3D();
		protected var wanderGlobal:SteerVector3D = new SteerVector3D();
		
		public var wanderStrength:Number = .1;
		public var wanderRate:Number = 0.8;
		public var wanderDirection:SteerVector3D = new SteerVector3D();
		
		public var slowingDistance:Number = 7;
		
		//
		protected var s1:SteerVector3D = new SteerVector3D();
		protected var s2:SteerVector3D = new SteerVector3D();
		protected var local:SteerVector3D = new SteerVector3D();
		protected var followTarget:SteerVector3D = new SteerVector3D();
		
		public var leaderAvoidWidth:Number = 4 * 0.5;
		public var leaderAvoidLength:Number = 30 * 0.5;
		public var fastslow:Boolean;
		
		public var leaderTarget:SteerVector3D = new SteerVector3D();
		public var leaderTargetCenter:SteerVector3D = new SteerVector3D();
		public var newTargetOffset:SteerVector3D = new SteerVector3D();
		
		public var leaderTargetRange:Number = 8;
		public var leaderTargetCount:Number = 0;
		
		public var vehicleList:Vehicle;
		public var vehicleLeader:Vehicle;
		
		public function LeaderFollowing(vehicleList:Vehicle = null, vehicleLeader:Vehicle = null)
		{
			this.vehicleList = vehicleList;
			this.vehicleLeader = vehicleLeader;
			
			fastslow = true;
		}

		override public function apply(veh:Vehicle):void
		{
			if(veh == vehicleLeader)
			{
				steeringForLeader(veh, accumulator);
			} else
			{
				steeringForFollower(veh, vehicleLeader, accumulator);
			}
			
			veh.applyGlobalForce(accumulator);
		}
		
		public function steeringForFollower(veh:Vehicle, target:Vehicle, v:SteerVector3D):void
		{
			steeringToFollowLeaderOrGetOutOfItsWay(veh, target, s1);
			steeringForSeparation(veh, 8 * veh.vehicleRadius, s2);
			
			v.x = s2.x * 0.2 + s1.x;
			v.y = s2.y * 0.2 + s1.y;
			v.z = s2.z * 0.2 + s1.z;
		}
		
		public function steeringToFollowLeaderOrGetOutOfItsWay(veh:Vehicle, target:Vehicle, v:SteerVector3D):void
		{
			target.localizePosition(veh.position, local);
			var adjust:Number = target.velocity.approximateLength() / target.maxSpeed;
			var margin:Number = 3 * veh.vehicleRadius;
			
			if(local.z > 0.0 && local.z < leaderAvoidLength * adjust + margin * 2.0 && local.x < leaderAvoidWidth * adjust + margin && local.x > -leaderAvoidWidth * adjust - margin && local.y < leaderAvoidWidth * adjust + margin && local.y > -leaderAvoidWidth * adjust - margin)
			{
				local.x += local.x <= 0 ? -1 : 1;
				local.y += local.y <= 0 ? -1 : 1;
				local.x *= 3;
				local.y *= 3;
				local.z *= 0.5;
				target.globalizePosition(local, followTarget);
				steeringForSeek(veh, followTarget, v);
			} else
			{
				var followDistance:Number = veh.vehicleRadius * -5;
				followTarget.x = target.forward.x * followDistance + target.position.x;
				followTarget.y = target.forward.y * followDistance + target.position.y;
				followTarget.z = target.forward.z * followDistance + target.position.z;
				
				steeringForArrival(veh, followTarget, v);
			}
		}
		
		public function steeringForLeader(veh:Vehicle, v:SteerVector3D):void
		{
			maybeSteerTowardRandomLeaderTarget(veh, v);
			if(v.equalsZero())
			{
				if(SteerVector3D.random() < 0.8)
				{
					var r:Number = SteerVector3D.random();
					if(r > 0.97) fastslow = true;
					if(r < 0.03) fastslow = false;
					var speed:Number = (fastslow ? 0.3 : 1) * veh.maxSpeed;
					accelToTargetSpeed(veh, speed, v);
				} 
				else {
					steeringForWander(veh, v);
				}
			}
		}
		
		public function maybeSteerTowardRandomLeaderTarget(veh:Vehicle, v:SteerVector3D):void
		{
			if(leaderTargetCount == 0 || leaderTarget.approximateDistance(veh.position) < veh.vehicleRadius * 3)
			{
				var loop:int = -1;
				var fwd:SteerVector3D = veh.forward;
				var px:Number = veh.position.x, py:Number = veh.position.y, pz:Number = veh.position.z;
				do
				{
					leaderTarget.x = leaderTargetRange * (SteerVector3D.random() * 0.8 + 0.1);
					leaderTarget.y = leaderTargetRange * (SteerVector3D.random() * 0.8 + 0.1);
					leaderTarget.z = leaderTargetRange * (SteerVector3D.random() * 0.8 + 0.1);
					
					newTargetOffset.x = leaderTarget.x - px;
					newTargetOffset.y = leaderTarget.y - py;
					newTargetOffset.z = leaderTarget.z - pz;
					
					newTargetOffset.setApproximateNormalize();
				} 
				while(newTargetOffset.dot(fwd) < -0.8 && ++loop < 50);
				
				leaderTargetCount = 50;
			} else
			{
				leaderTargetCount--;
			}
			
			if(SteerVector3D.random() < 0.43)
			{
				steeringForSeek(veh, leaderTarget, v);
			} else
			{
				v.setZero();
			}
		}
		
		public function accelToTargetSpeed(veh:Vehicle, targetSpeed:Number, v:SteerVector3D):void
		{
			var speed:Number = veh.velocity.approximateLength();
			if(speed < targetSpeed) 
			{
				v.setScale(veh.maxForce, veh.forward);
			} else {
				v.setScale(-veh.maxForce, veh.forward);
			}
		}

		public function steeringForSeek(veh:Vehicle, target:SteerVector3D, v:SteerVector3D):void
		{
			v.setDiff(target, veh.position);
			v.setApproximateTruncate(veh.maxSpeed);
			v.setDiff(v, veh.velocity);
		}
		
		public function steeringForWander(veh:Vehicle, v:SteerVector3D):void
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
			
			v.x = veh.forward.x * SQRT2 + wanderGlobal.x;
			v.y = veh.forward.y * SQRT2 + wanderGlobal.y;
			v.z = veh.forward.z * SQRT2 + wanderGlobal.z;
		}
		
		public function steeringForSeparation(veh:Vehicle, separationDistance:Number, v:SteerVector3D):void
		{			
			var other:Vehicle = vehicleList;
			var approximateDistance:Number;
			
			var px:Number = veh.position.x;
			var py:Number = veh.position.y;
			var pz:Number = veh.position.z;
			var dx:Number, dy:Number, dz:Number;
			var vx:Number = 0, vy:Number = 0, vz:Number = 0;
			var a:Number, b:Number, c:Number, t:Number;
			var f:Number;
			
			while(other) 
			{
				if(other == veh)
				{
					other = other.next;
					continue;
				}
				dx = px - other.position.x;
				dy = py - other.position.y;
				dz = pz - other.position.z;
				
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
				
				approximateDistance = a * 0.9375 + (b + c) * 0.375;
				
				if( approximateDistance < separationDistance && approximateDistance > 0) 
				{					
					
					f = 1 / (approximateDistance * approximateDistance);
					
					vx += dx * f;
					vy += dy * f;
					vz += dz * f;
				}
				other = other.next;
			}
			
			v.x = vx;
			v.y = vy;
			v.z = vz;
		}
		
		public function steeringForArrival(veh:Vehicle, target:SteerVector3D, v:SteerVector3D):void
		{
			var dx:Number = target.x - veh.position.x;
			var dy:Number = target.y - veh.position.y;
			var dz:Number = target.z - veh.position.z;
			
			var a:Number = dx;
			if (a < 0) a = -a;
			
			var b:Number = dy;
			if (b < 0) b = -b;
			
			var c:Number = dz;
			if (c < 0) c = -c;
			
			var t:Number;
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
			
			var distance:Number =  a * 0.9375 + (b + c) * 0.375;
			
			var rampedSpeed:Number = veh.maxSpeed * (distance / slowingDistance);
			var clippedSpeed:Number = (rampedSpeed > veh.maxSpeed ? veh.maxSpeed : rampedSpeed) / distance;
			
			v.x = dx * clippedSpeed - veh.velocity.x;
			v.y = dy * clippedSpeed - veh.velocity.y;
			v.z = dz * clippedSpeed - veh.velocity.z;
		}
	}
}
