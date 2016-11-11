package ru.inspirit.steering 
{
	import ru.inspirit.steering.behavior.AbstractBehavior;

	/**
	 * Vehicle Group
	 * simple way to control large number of vehicles
	 * 
	 * @author Eugene Zatepyakin
	 */
	public class VehicleGroup 
	{
		public var vehicleFirst:Vehicle;
		public var vehicleLast:Vehicle;
		
		public var defaultBehavior:AbstractBehavior;
		
		/**
		 * @param defaultBehavior	is used if Vehicle doesn't have any
		 */
		public function VehicleGroup(defaultBehavior:AbstractBehavior)
		{
			this.defaultBehavior = defaultBehavior;
		}
		
		/**
		 * Update all vehicles in group
		 */
		public function update():void
		{
			var veh:Vehicle = vehicleFirst;
			var behave:AbstractBehavior;
			
			while(veh) 
			{
				behave = veh.behaviorList ? veh.behaviorList.behaviorFirst : defaultBehavior;
				while(behave)
				{
					behave.apply(veh);
					behave = behave.next;
				}
				
				veh.update();
				
				veh = veh.next;
			}
		}
		
		/**
		 * Add vehicle to group
		 * 
		 * @param veh	Vehicle object to add
		 */
		public function addVehicle(veh:Vehicle):void
		{
			if(vehicleFirst)
			{
				vehicleLast = vehicleLast.next = veh;
			}
			else 
			{
				vehicleFirst = vehicleLast = veh;
			}
		}
		
		/**
		 * Remove vehicle from group
		 * 
		 * @param veh	Vehicle to remove
		 */
		public function removeVehicle(veh:Vehicle):void
		{
			var node:Vehicle = vehicleFirst;
            if (veh == vehicleFirst)
            {                    
                if (vehicleFirst == vehicleLast)
                {
					vehicleFirst = vehicleLast = null;
                }
                else
                {
                    vehicleFirst = vehicleFirst.next;
                    node.next = null;
                    if (vehicleFirst == null) vehicleLast = null;
                }
            } 
            else 
            {
	            while (node.next != veh) node = node.next;
	            if (node.next == vehicleLast) vehicleLast = node;
	            node.next = veh.next;
			}
		}
		
		/**
		 * Remove all Vehicles
		 */
		public function clear():void
        {
            var node:Vehicle = vehicleFirst;
            vehicleFirst = null;
            
            var next:Vehicle;
            while (node)
            {
                next = node.next;
                node.next = null;
                node = next;
            }
        }
	}
}
