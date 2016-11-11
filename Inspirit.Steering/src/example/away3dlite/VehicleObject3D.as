package example.away3dlite 
{
	import flash.geom.Vector3D;
	import away3dlite.core.base.Object3D;
	import ru.inspirit.steering.SteerVector3D;
	import ru.inspirit.steering.Vehicle;
	import ru.inspirit.steering.behavior.BehaviorList;

	/**
	 * @author Eugene Zatepyakin
	 */
	public class VehicleObject3D extends Vehicle 
	{
		public static var SCALE:Number = 10;
		
		public var obj3d:Object3D;
		public var lookTarget:Vector3D = new Vector3D();
		
		public function VehicleObject3D(obj3d:Object3D, position:SteerVector3D = null, behaviorList:BehaviorList = null)
		{
			super(position, behaviorList);
			
			this.obj3d = obj3d;
		}
		
		override public function update():void
		{
			super.update();
			
			lookTarget.x = obj3d.x = position.x * SCALE;
			lookTarget.y = obj3d.y = position.y * SCALE;
			lookTarget.z = obj3d.z = position.z * SCALE;
			
			lookTarget.incrementBy(velocity);
			
			obj3d.lookAt(lookTarget, Vector3D.Y_AXIS);
		}
	}
}
