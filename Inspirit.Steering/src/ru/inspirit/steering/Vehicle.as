package ru.inspirit.steering 
{
	import ru.inspirit.steering.behavior.BehaviorList;

	import flash.geom.Vector3D;

	/**
	 * Base class for any moving entity
	 * 
	 * @author Eugene Zatepyakin
	 */
	public class Vehicle 
	{		
		public static const globalUp:SteerVector3D = new SteerVector3D(0, 0.1, 0);
		
		// BOUNDS OPTIONS
		public static const EDGE_NONE:uint = 0;
		public static const EDGE_WRAP:uint = 1;
		public static const EDGE_BOUNCE:uint = 2;
		
		public var boundsCentre:SteerVector3D = new SteerVector3D();
		public var boundsRadius:Number = 25;
		public var edgeBehavior:uint = EDGE_WRAP;
		
		public var mass:Number;
		public var maxSpeed:Number;
		public var maxForce:Number;
		public var accelDamping:Number;
		
		public var velocity:SteerVector3D;
		public var allForces:SteerVector3D;
		public var acceleration:SteerVector3D;
		
		public var position:SteerVector3D;
		public var forward:SteerVector3D;
		public var side:SteerVector3D;
		public var up:SteerVector3D;
		
		public var vehicleRadius:Number;
		
		public var next:Vehicle;
		public var behaviorList:BehaviorList;
		
		public function Vehicle(position:SteerVector3D = null, behaviorList:BehaviorList = null) 
		{
			identity();
			if (position) {
				this.position = position;
			}
			
			this.behaviorList = behaviorList;
		}

		public function identity():void
		{
			mass = 1;
			maxSpeed = 1;
			maxForce = 0.04;
			accelDamping = 0.7;
			vehicleRadius = 3;
			
			acceleration = new SteerVector3D();
			velocity = new SteerVector3D();
			allForces = new SteerVector3D();
			
			position = new SteerVector3D();
			forward = new SteerVector3D(0, 0, 1);
			side = new SteerVector3D(1, 0, 0);
			up = new SteerVector3D(0, 1, 0);
		}
		
		public function applyGlobalForce(force:SteerVector3D):void
		{
			allForces.x += force.x;
			allForces.y += force.y;
			allForces.z += force.z;
		}

		public function update():void
		{
			allForces.setApproximateTruncate(maxForce);
			
			var nax:Number, nay:Number, naz:Number;
			var m:Number;
			
			if(mass == 1) 
			{
				nax = allForces.x;
				nay = allForces.y;
				naz = allForces.z;
			} else {
				m = 1 / mass;
				nax = allForces.x * m;
				nay = allForces.y * m;
				naz = allForces.z * m;
			}
			
			allForces.x = allForces.y = allForces.z = 0;
			
			acceleration.x = nax + accelDamping * (acceleration.x - nax);
			acceleration.y = nay + accelDamping * (acceleration.y - nay);
			acceleration.z = naz + accelDamping * (acceleration.z - naz);
			
			velocity.x += acceleration.x;
			velocity.y += acceleration.y;
			velocity.z += acceleration.z;
			
			velocity.setApproximateTruncate(maxSpeed);
			
			position.x += velocity.x;
			position.y += velocity.y;
			position.z += velocity.z;
			
			boundsCheck();
			
			var aux:Number = acceleration.x * 0.5 + up.x + globalUp.x;
			var auy:Number = acceleration.y * 0.5 + up.y + globalUp.y;
			var auz:Number = acceleration.z * 0.5 + up.z + globalUp.z;
			
			m = Math.sqrt( aux * aux + auy * auy + auz * auz );
			if(m != 0) 
			{
				m = 1 / m;
				aux *= m;
				auy *= m;
				auz *= m;
			}
			
			var speed:Number = velocity.magnitude();
			if(speed > 0) 
			{
				speed = 1 / speed;
				forward.x = velocity.x * speed;
				forward.y = velocity.y * speed;
				forward.z = velocity.z * speed;
				
				side.x = forward.y * auz - forward.z * auy;
				side.y = forward.z * aux - forward.x * auz;
				side.z = forward.x * auy - forward.y * aux;
				
				up.setCross(side, forward);
			}
		}
		
		public function boundsCheck():void
		{
			if(edgeBehavior == EDGE_NONE) return;
			
			var distance:Number = Vector3D.distance(position, boundsCentre);
				
			if( distance > boundsRadius + vehicleRadius ) 
			{
				switch( edgeBehavior )
				{
					case EDGE_BOUNCE :
						
						position.decrementBy(boundsCentre);
						position.normalize();
						position.scaleBy(boundsRadius + vehicleRadius);
						
						velocity.scaleBy(-1);
						position.incrementBy(velocity);
						position.incrementBy(boundsCentre);
					
						break;
					
					case EDGE_WRAP :

						//position.decrementBy(boundsCentre);
						//position.negate();
						//position.incrementBy(boundsCentre);
						var wrap:SteerVector3D = position.scale(1 / distance);
						wrap.scaleBy(-2 * boundsRadius);
						position.incrementBy(wrap);

						break;
				}
			}
		}

		public function globalizePosition(local:SteerVector3D, globalized:SteerVector3D):void
		{
			globalized.x = (side.x * local.x) + (up.x * local.y) + (forward.x * local.z) + position.x;
			globalized.y = (side.y * local.x) + (up.y * local.y) + (forward.y * local.z) + position.y;
			globalized.z = (side.z * local.x) + (up.z * local.y) + (forward.z * local.z) + position.z;
		}
		
		public function globalizeDirection(local:SteerVector3D, globalized:SteerVector3D):void
		{	
			globalized.x = (side.x * local.x) + (up.x * local.y) + (forward.x * local.z);
			globalized.y = (side.y * local.x) + (up.y * local.y) + (forward.y * local.z);
			globalized.z = (side.z * local.x) + (up.z * local.y) + (forward.z * local.z);
		}
		
		public function localizePosition(global:SteerVector3D, localized:SteerVector3D):void
		{
			//component.setDiff(global, position);
			var xd:Number = global.x - position.x;
			var yd:Number = global.y - position.y;
			var zd:Number = global.z - position.z;
			//localized.setTo(component.dot(side), component.dot(up), component.dot(forward));
			localized.x = xd * side.x + yd * side.y + zd * side.z;
			localized.y = xd * up.x + yd * up.y + zd * up.z;
			localized.z = xd * forward.x + yd * forward.y + zd * forward.z;
		}
	}
}
