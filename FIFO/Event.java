
// event representation
class Event implements Comparable {

 public Event(int a_type, double a_time, int order) { _type = a_type; time = a_time; order_num = order;}
  
 public double time;
 private int order_num;
 private int _type;
 public double arrive_router_time;
 
 public int get_order() { return order_num; }
 public int get_type() { return _type; }
 public double get_time() { return time; }
 public void set_time(double a_time) {arrive_router_time = a_time;}
 public double get_arrive_router_time() {return arrive_router_time;}
 public Event leftlink, rightlink, uplink;
 
 
 public int compareTo(Object _cmpEvent ) {
  double _cmp_time = ((Event) _cmpEvent).get_arrive_router_time() ;
  if( this.arrive_router_time < _cmp_time) return -1;
  if( this.arrive_router_time == _cmp_time) return 0;
  return 1;
 }
};
