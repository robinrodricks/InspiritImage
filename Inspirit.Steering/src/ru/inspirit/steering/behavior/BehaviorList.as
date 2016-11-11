package ru.inspirit.steering.behavior 
{

	/**
	 * Behavior list is based on Linked List structure
	 * Every Vehicle have a list of behaviors to apply
	 * 
	 * @author Eugene Zatepyakin
	 */
	public final class BehaviorList 
	{
		public var behaviorFirst:AbstractBehavior;
		public var behaviorLast:AbstractBehavior;
		
		public function BehaviorList(...args)
		{
			if(args.length > 0) addBehavior.apply(this, args);
		}
		
		/**
		 * @param args	comma separated list of Behaviors to add
		 */
		public function addBehavior(...args):void
		{
			var k:int = args.length;
			
			var node:AbstractBehavior = args[0];
			
			if(behaviorFirst)
			{
				behaviorLast = behaviorLast.next = node;
			}
			else 
			{
				behaviorFirst = behaviorLast = node;
			}
			
			if(k > 1)
			{
				for(var i:int = 1; i < k; ++i)
				{
					node = args[i];
					behaviorLast = behaviorLast.next = node;
				}
			}
		}
		
		public function removeBehavior(behave:AbstractBehavior):void
		{
			var node:AbstractBehavior = behaviorFirst;
            if (behave == behaviorFirst)
            {                    
                if (behaviorFirst == behaviorLast)
                {
					behaviorFirst = behaviorLast = null;
                }
                else
                {
                    behaviorFirst = behaviorFirst.next;
                    node.next = null;
                    if (behaviorFirst == null) behaviorLast = null;
                }
            } 
            else 
            {
	            while (node.next != behave) node = node.next;
	            if (node.next == behaviorLast) behaviorLast = node;
	            node.next = behave.next;
			}
		}
		
		public function clear():void
        {
            var node:AbstractBehavior = behaviorFirst;
            behaviorFirst = null;
            
            var next:AbstractBehavior;
            while (node)
            {
                next = node.next;
                node.next = null;
                node = next;
            }
        }
	}
}
