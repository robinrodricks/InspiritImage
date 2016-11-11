package example.ribbon 
{
	import ru.inspirit.steering.Vehicle;
	import ru.inspirit.steering.VehicleGroup;
	import ru.inspirit.steering.behavior.AbstractBehavior;
	import ru.inspirit.steering.behavior.BehaviorList;
	import ru.inspirit.steering.behavior.Wander;
	import ru.inspirit.steering.behavior.combined.Flocking;
	import ru.inspirit.steering.behavior.combined.LeaderFollowing;

	import flash.display.Shape;
	import flash.display.Sprite;
	import flash.events.Event;
	import flash.geom.Rectangle;
	import flash.utils.getTimer;

	/**
	 * @author Eugene Zatepyakin
	 */
	public class ExampleRibbon extends Sprite 
	{
		public static const SCALE:uint = 10;
		
		public var numberOfVehicles:uint = 10;
		
		public var viewport:Shape;
		public var ribbManager:Ribbon3DManager;
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
		
		public function ExampleRibbon()
		{
			if(stage) init();
			else addEventListener(Event.ADDED_TO_STAGE, init);
		}
		
		private function init(e:Event = null):void
		{
			removeEventListener(Event.ADDED_TO_STAGE, init);
			
			viewport = new Shape();
			viewport.x = 400;
			viewport.y = 300;
			
			ribbManager = new Ribbon3DManager(viewport);
			
			addChild(viewport);
			
			scrollRect = new Rectangle(0, 0, 800, 600);
		}
		
		public function initFlockingBoids(num:uint = 10):void
		{
			numberOfVehicles = num;
			
			var flockBehave:Flocking = new Flocking();
			
			createVehicles(flockBehave, Vehicle.EDGE_BOUNCE);
			
			flockBehave.separate.separateList = flockBehave.align.alignList = flockBehave.cohere.cohereList = vehGroup.vehicleFirst;
		}
		
		public function initFollowers(num:uint = 10):void
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
			var veh:Ribbon3D;
			var n:int = numberOfVehicles;
			while(--n > -1)
			{				
				veh = new Ribbon3D(ribbManager, null, new BehaviorList(behave), int(Math.random()*0xFFFFFF), int(Math.random()*3)+2, 100+int(Math.random()*100));
				veh.vehicleRadius = 0.5;
				veh.maxSpeed = 0.6;
				veh.maxForce = 0.07;
				veh.edgeBehavior = edgeBehavior;
				veh.position.setUnitRandom();
				veh.position.scaleBy(SCALE);
				veh.velocity.setUnitRandom();
				veh.velocity.scaleBy(veh.maxSpeed);
				vehGroup.addVehicle(veh);
			}
		}
		
		public function render(e:Event = null):void
		{
			var t:uint = getTimer();
			
			ribbManager.prerender();
			//viewport.graphics.clear();
			
			vehGroup.update();
			
			steer_procc_t += getTimer() - t;
			steer_procc_n ++;
			
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
