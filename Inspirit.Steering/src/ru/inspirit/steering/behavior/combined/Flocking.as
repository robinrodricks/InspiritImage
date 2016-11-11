package ru.inspirit.steering.behavior.combined 
{
	import ru.inspirit.steering.Vehicle;
	import ru.inspirit.steering.behavior.AbstractBehavior;
	import ru.inspirit.steering.behavior.Alignment;
	import ru.inspirit.steering.behavior.Cohesion;
	import ru.inspirit.steering.behavior.Separation;

	/**
	 * Flocking Behavior
	 * as described by Craig Reynolds 
	 * 
	 * @author Eugene Zatepyakin
	 */
	public final class Flocking extends AbstractBehavior 
	{
		public var separate:Separation;
		public var align:Alignment;
		public var cohere:Cohesion;
		
		public function Flocking(separateList:Vehicle = null, alignList:Vehicle = null, cohereList:Vehicle = null)
		{
			separate = new Separation(separateList);
			align = new Alignment(alignList);
			cohere = new Cohesion(cohereList);
		}

		override public function apply(veh:Vehicle):void
		{
			separate.apply(veh);
			align.apply(veh);
			cohere.apply(veh);
		}
	}
}
