package example.away3dlite 
{
	import away3dlite.cameras.Camera3D;
	import away3dlite.containers.Scene3D;
	import away3dlite.containers.View3D;
	import away3dlite.core.clip.RectangleClipping;
	import away3dlite.core.render.BasicRenderer;
	import away3dlite.materials.ColorMaterial;
	import away3dlite.primitives.Cone;

	import ru.inspirit.steering.Vehicle;
	import ru.inspirit.steering.VehicleGroup;
	import ru.inspirit.steering.behavior.AbstractBehavior;
	import ru.inspirit.steering.behavior.BehaviorList;
	import ru.inspirit.steering.behavior.Wander;
	import ru.inspirit.steering.behavior.combined.Flocking;
	import ru.inspirit.steering.behavior.combined.LeaderFollowing;

	import flash.display.Sprite;
	import flash.events.Event;
	import flash.utils.getTimer;

	/**
	 * @author Eugene Zatepyakin
	 */
	public class ExampleFlocking extends Sprite 
	{
		public static const SCALE:uint = 10;
		
		public var scene:Scene3D;
		public var camera:Camera3D;
		public var view:View3D;
		public var renderer:BasicRenderer = new BasicRenderer();
		public var clipping:RectangleClipping = new RectangleClipping();
		
		public var numberOfVehicles:uint = 150;
		
		public var cones:Vector.<Cone> = new Vector.<Cone>(numberOfVehicles, true);
		public var vehGroup:VehicleGroup;
		
		public var steerTime:int = 0;
		public var renderTime:int = 0;
		
		protected var steer_procc_t:uint = 0;
		protected var render_procc_t:uint = 0;
		protected var steer_procc_n:uint = 0;
		protected var render_procc_n:uint = 0;
		
		protected var _timer:uint;
		protected var _fps:uint;
		protected var _ms:uint;
		protected var _ms_prev:uint;

		public function ExampleFlocking()
		{
			if(stage) init();
			else addEventListener(Event.ADDED_TO_STAGE, init);
		}
		
		private function init(e:Event = null):void
		{
			removeEventListener(Event.ADDED_TO_STAGE, init);
			
			scene = new Scene3D();
			
			camera = new Camera3D();
			camera.z = -800;
			
			view = new View3D();
			view.scene = scene;
			view.camera = camera;
			view.renderer = renderer;
			view.clipping = clipping;
			view.x = 400;
			view.y = 300;
			
			addChild(view);
		}
		
		public function initFlockingBoids(num:uint = 100):void
		{
			numberOfVehicles = num;
			
			var flockBehave:Flocking = new Flocking();
			
			createVehicles(flockBehave, Vehicle.EDGE_WRAP);
			
			flockBehave.separate.separateList = flockBehave.align.alignList = flockBehave.cohere.cohereList = vehGroup.vehicleFirst;
		}
		
		public function initFollowers(num:uint = 100):void
		{
			numberOfVehicles = num;
			
			var followBehave:LeaderFollowing = new LeaderFollowing();
			
			createVehicles(followBehave);
			
			followBehave.vehicleList = vehGroup.vehicleFirst;
			followBehave.vehicleLeader = vehGroup.vehicleFirst;
			followBehave.leaderAvoidWidth = 4 * vehGroup.vehicleFirst.vehicleRadius;
			followBehave.leaderAvoidLength = 30 * vehGroup.vehicleFirst.vehicleRadius;
		}
		
		public function createVehicles(behave:AbstractBehavior, edgeBehavior:uint = Vehicle.EDGE_NONE):void
		{
			vehGroup = new VehicleGroup(new Wander());
			var veh:VehicleObject3D;
			var cone:Cone;
			var n:int = numberOfVehicles;
			var i:int = -1;
			while(--n > -1)
			{
				cone = new Cone(new ColorMaterial(int(Math.random()*0xFFFFFF)), 7, 14, 3, 1, false, false);
				cone.bothsides = true;
				cones[int(++i)] = cone;
				
				veh = new VehicleObject3D(cone, null, new BehaviorList(behave));
				veh.vehicleRadius = 1;
				veh.maxSpeed = 0.6;
				veh.maxForce = 0.07;
				veh.edgeBehavior = edgeBehavior;
				veh.position.setUnitRandom();
				veh.position.scaleBy(SCALE);
				veh.velocity.setUnitRandom();
				veh.velocity.scaleBy(veh.maxSpeed);
				vehGroup.addVehicle(veh);
				
				scene.addChild(cone);
			}
		}
		
		public function render(e:Event = null):void
		{
			var t:uint = getTimer();
			
			vehGroup.update();
			
			steer_procc_t += getTimer() - t;
			steer_procc_n ++;
			
			t = getTimer();
			
			view.render();
			
			render_procc_t += getTimer() - t;
			render_procc_n ++;
			
			countFrameTime();
		}
		
		public function countFrameTime():void
		{
			_timer = getTimer();
			if( _timer - 1000 >= _ms_prev )
			{
				_ms_prev = _timer;

				steerTime = int(steer_procc_t / steer_procc_n + 0.5);
				renderTime = int(render_procc_t / render_procc_n + 0.5);

				steer_procc_n = steer_procc_t = render_procc_n = render_procc_t = _fps = 0;
			}

			_fps ++;
			_ms = _timer;
		}
	}
}
