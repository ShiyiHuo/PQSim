import java.util.*;
class Sim {

// Class Sim variables
public static double Clock, MeanInterArrivalTime, MeanServiceTime, SourceServiceTime, RouterDelayTime, SourceDelayTime, SIGMA, LastEventTime,
        TotalBusy, MaxQueueLength, SumResponseTime, NormalMean, NormalStandarDeviation, TotalDelay;
public static long  NumberOfCustomers, QueueLength, NumberInService,
        TotalCustomers, NumberOfDepartures, LongService;
public static int RouterArrivalRate, RouterServiceRate, SourceServiceRate;

public final static int arrival = 1;
public final static int departure = 2;
public static ArrayList<Event> packages;
public static ArrayList<Event> finalArrivals;
public static EventList FutureEventList;
public static Queue Customers;
public static Random stream;
public static Event finished;
public static int dropped;
public static double totalWait = 0;
public static int count;

public static void main(String argv[]) {
  TotalDelay = 0;
 // RouterArrivalRate = Double.parseDouble(argv[1]);
  RouterServiceRate = 1250;
  SourceServiceRate = 1250;
  SourceServiceTime = 1/SourceServiceRate;
  NormalMean = Double.parseDouble(argv[1]);
  NormalStandarDeviation = Double.parseDouble(argv[2]);
  
  MeanInterArrivalTime = 1.0/1125; MeanServiceTime = 1/RouterServiceRate;
  SIGMA                = 0.6; TotalCustomers  = 300000;
  long seed            = Long.parseLong(argv[0]);

  stream = new Random(seed);           // initialize rng stream
  RouterDelayTime = 0.05;
  SourceDelayTime = normal(stream, NormalMean, NormalStandarDeviation);
  FutureEventList = new EventList();
  packages = new ArrayList<Event>();
  Customers = new Queue();
  count = 0;
  Initialization();

  // Loop until first "TotalCustomers" have departed
  while(NumberOfDepartures < TotalCustomers ) {
    Event evt = (Event)FutureEventList.getMin();  // get imminent event
    FutureEventList.dequeue();                    // be rid of it
    Clock = evt.get_time();                       // advance simulation time
    if( evt.get_type() == arrival ) ProcessArrival(evt);
    else  ProcessDeparture(evt);
    }
	System.out.println("size "+packages.size());

	Collections.sort(packages);
	int unOrderedCount = 0;
	
	int currSequence = -1;
	
	for(int i=0; i<packages.size(); i++){
	//	System.out.println(packages.get(i).get_arrive_router_time());
		if(i==0){
			currSequence = packages.get(i).get_order();
			continue;
		}
		if(packages.get(i).get_order()>currSequence){
			currSequence = packages.get(i).get_order();
		} else{
			unOrderedCount++;
		}
	}
//	System.out.println("maxqueuelength "+MaxQueueLength);
//	System.out.println("the ratio of out of order packets is "+unOrderedCount/200000.0);
//	System.out.println("the average package delay is "+(TotalDelay/200000.0+0.05));
    Clock = packages.get(0).get_arrive_router_time();
	System.out.println("Clock is : "+Clock);
	LastEventTime = Clock;
	count = 0;
	NumberInService = 0;
	currSequence = -1;
	NumberOfDepartures = 0;
	FutureEventList = new EventList();
	Event e = new Event(arrival, packages.get(0).get_arrive_router_time(), count);
//	HighQueue.enqueue(e);
	FutureEventList.enqueue(e);
	int j = 0;
	ArrayList<Event> arr = new ArrayList<Event>();
	finalArrivals = new ArrayList<Event>();
	Customers = new Queue();
	dropped = 0;
	QueueLength = 0;
	MaxQueueLength = 0;
  // Start to process the arrivals and departures on router
  while(NumberOfDepartures < TotalCustomers ) {
	//System.out.println("");
    Event evt = (Event)FutureEventList.getMin();  // get imminent event
    FutureEventList.dequeue();                    // be rid of it
    Clock = evt.get_time();                       // advance simulation time
    
	if(j==0){
		
		ProcessRouterArrival(evt);
		j++;
		continue;
	}
    if( evt.get_type() == arrival ) {
	//	System.out.println("Arrival at time: "+evt.get_time());
		if(QueueLength>=10000){
			dropped++;
			count++;
		  if(count<TotalCustomers){
			  Event next_arrival = new Event(arrival, packages.get(count).get_arrive_router_time(), packages.get(count).get_order());
			  FutureEventList.enqueue( next_arrival );
			  
			  LastEventTime = Clock; 
		  }
		}else{
			ProcessRouterArrival(evt);
		}
		
	}
    else  {
	//	System.out.println("Departure at time: "+evt.get_time());
		ProcessRouterDeparture(evt);
	}		
	j++;
    }

	unOrderedCount = 0;
	//Collections.sort(finalArrivals);

	for(int i=0; i<finalArrivals.size(); i++){
		
		if(i==0){
			currSequence = finalArrivals.get(i).get_order();
			continue;
		}
		if(finalArrivals.get(i).get_order()>currSequence){
			currSequence = finalArrivals.get(i).get_order();
		} else{
			unOrderedCount++;
		}
	}
	System.out.println("the ratio of out of order packets is "+unOrderedCount/(double)TotalCustomers);
	System.out.println("the average package delay is "+(TotalDelay/(double)TotalCustomers));
	System.out.println("the ratio of dropped packets is "+(dropped/(double)TotalCustomers));
	System.out.println("max queue length is "+MaxQueueLength);
	System.out.println("Avg wait time is "+totalWait/(double)TotalCustomers);

  ReportGeneration();
 }

 // seed the event list with TotalCustomers arrivals
 public static void Initialization()   { 
  Clock = 0.0;
  QueueLength = 0;
  NumberInService = 0;
  LastEventTime = 0.0;
  TotalBusy = 0 ;
  MaxQueueLength = 0;
  SumResponseTime = 0;
  NumberOfDepartures = 0;
  LongService = 0;

  // create first arrival event
  Event evt = new Event(arrival, exponential( stream, MeanInterArrivalTime), count);
  count++;
  FutureEventList.enqueue( evt );
 }
  public static void ProcessRouterArrival(Event evt) {

  Customers.enqueue(evt); 
  QueueLength++;
  // if the server is idle, fetch the event, do statistics
  // and put into service
  if( NumberInService == 0) {
	ScheduleRouterDeparture();
  }
  else TotalBusy += (Clock - LastEventTime);  // server is busy

  // adjust max queue length statistics
  if (MaxQueueLength < QueueLength) MaxQueueLength = QueueLength;

  // schedule the next arrival
//  System.out.println("count is "+count);
  count++;
  if(count<TotalCustomers){
	  Event next_arrival = new Event(arrival, packages.get(count).get_arrive_router_time(), packages.get(count).get_order());
	  FutureEventList.enqueue( next_arrival );
	  
	  LastEventTime = Clock; 
  }
 
 }

  public static void ScheduleRouterDeparture() {
  double ServiceTime = MeanServiceTime;
  // get the job at the head of the queue
 // while (( ServiceTime = normal(stream, MeanServiceTime, SIGMA)) < 0 );
  
  while (( ServiceTime = 1.0/1250.0) < 0 );
  Event evt = (Event)Customers.dequeue();
  QueueLength--;
  Event depart = new Event(departure,Clock+ServiceTime, evt.get_order());
  finished = depart;
  TotalDelay+=(finished.get_time()-evt.get_time());
  
  totalWait+=(finished.get_time()-evt.get_time()-1.0/1250.0);
  
  finished.set_time(evt.get_arrive_router_time());
  FutureEventList.enqueue( depart );
  
  NumberInService = 1;
 }
 

public static void ProcessRouterDeparture(Event e) {

 // get the customer description
 //Event finished = (Event) Customers.dequeue();
 double delay = 0.05;
 TotalDelay += delay;
 
 finalArrivals.add(finished);
 // if there are customers in the queue then schedule
 // the departure of the next one
  if( QueueLength > 0 ) ScheduleRouterDeparture();
  else NumberInService = 0;
  // measure the response time and add to the sum
  double response = (Clock - finished.get_time());
  SumResponseTime += response;
  if( response > 4.0 ) LongService++; // record long service
  TotalBusy += (Clock - LastEventTime );
  NumberOfDepartures++;
  LastEventTime = Clock;
 }

 public static void ProcessArrival(Event evt) {
  Customers.enqueue(evt); 
  QueueLength++;
  // if the server is idle, fetch the event, do statistics
  // and put into service
  if( NumberInService == 0) ScheduleDeparture();
  else TotalBusy += (Clock - LastEventTime);  // server is busy

  // adjust max queue length statistics
  if (MaxQueueLength < QueueLength) MaxQueueLength = QueueLength;

  // schedule the next arrival
  Event next_arrival = new Event(arrival, Clock+exponential(stream, MeanInterArrivalTime), count);
  count++;
  FutureEventList.enqueue( next_arrival );
  LastEventTime = Clock;
 }

 public static void ScheduleDeparture() {
  double ServiceTime = MeanServiceTime;
  // get the job at the head of the queue
  while (( ServiceTime = 1.0/1250.0) < 0 );
 // ServiceTime = SourceServiceTime;
  Event depart = new Event(departure,Clock+ServiceTime, count);
  FutureEventList.enqueue( depart );
  NumberInService = 1;
  QueueLength--;
 }

public static void ProcessDeparture(Event e) {
 // get the customer description
 Event finished = (Event) Customers.dequeue();
 double delay = normal(stream, NormalMean, NormalStandarDeviation);
 totalWait+=(Clock-finished.get_time()-1.0/1250.0);
 TotalDelay += (Clock+delay-finished.get_time());
 finished.set_time(Clock+delay);
 packages.add(finished);
 // if there are customers in the queue then schedule
 // the departure of the next one
  if( QueueLength > 0 ) ScheduleDeparture();
  else NumberInService = 0;
  // measure the response time and add to the sum
  double response = (Clock - finished.get_time());
  SumResponseTime += response;
  if( response > 4.0 ) LongService++; // record long service
  TotalBusy += (Clock - LastEventTime );
  NumberOfDepartures++;
  LastEventTime = Clock;
 }

public static void ReportGeneration() {
double RHO   = TotalBusy/Clock;
double AVGR  = SumResponseTime/TotalCustomers;
double PC4   = ((double)LongService)/TotalCustomers;
/*

System.out.println( "SINGLE SERVER QUEUE SIMULATION - GROCERY STORE CHECKOUT COUNTER ");
System.out.println( "\tMEAN INTERARRIVAL TIME                         " 
	+ MeanInterArrivalTime );
System.out.println( "\tMEAN SERVICE TIME                              " 
	+ MeanServiceTime );
System.out.println( "\tSTANDARD DEVIATION OF SERVICE TIMES            " + SIGMA );
System.out.println( "\tNUMBER OF CUSTOMERS SERVED                     " + TotalCustomers );
System.out.println(); 
System.out.println( "\tSERVER UTILIZATION                             " + RHO );
System.out.println( "\tMAXIMUM LINE LENGTH                            " + MaxQueueLength );
System.out.println( "\tAVERAGE RESPONSE TIME                          " + AVGR + "  MINUTES" );
System.out.println( "\tPROPORTION WHO SPEND FOUR "); 
System.out.println( "\t MINUTES OR MORE IN SYSTEM                     " + PC4 );
System.out.println( "\tSIMULATION RUNLENGTH                           " + Clock + " MINUTES" );
System.out.println( "\tNUMBER OF DEPARTURES                           " + TotalCustomers );*/
}

public static double exponential(Random rng, double mean) {
 return -mean*Math.log( rng.nextDouble() );
}

public static double SaveNormal;
public static int  NumNormals = 0;
public static final double  PI = 3.1415927 ;

public static double normal(Random rng, double mean, double sigma) {
        double ReturnNormal;
        // should we generate two normals?
        if(NumNormals == 0 ) {
          double r1 = rng.nextDouble();
          double r2 = rng.nextDouble();
          ReturnNormal = Math.sqrt(-2*Math.log(r1))*Math.cos(2*PI*r2);
          SaveNormal   = Math.sqrt(-2*Math.log(r1))*Math.sin(2*PI*r2);
          NumNormals = 1;
        } else {
          NumNormals = 0;
          ReturnNormal = SaveNormal;
        }
        return ReturnNormal*sigma + mean ;
 }
}

