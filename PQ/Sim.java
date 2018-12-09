import java.util.*;
class Sim {

// Class Sim variables
public static double Clock, MeanInterArrivalTime, MeanServiceTime, SourceServiceTime, RouterDelayTime, SourceDelayTime, SIGMA, LastEventTime,
        TotalBusy, MaxQueueLength, MaxQueueLength1, MaxQueueLength2, SumResponseTime, NormalMean, NormalStandarDeviation, TotalDelay;
public static long  NumberOfCustomers, QueueLength, LowQueueLength, HighQueueLength, NumberInService,
        TotalCustomers, NumberOfDepartures, LongService;
public static int RouterArrivalRate, RouterServiceRate, SourceServiceRate;

public final static int arrival = 1;
public final static int departure = 2;
public static ArrayList<Event> packages;
public static ArrayList<Event> finalArrivals;
public static EventList FutureEventList;
public static Queue HighQueue;
public static Queue LowQueue;
public static Queue Customers;
public static Random stream;
public static Event finished;

public static int totalWait = 0;
public static int count;

/******* newly added ******/
public static double[] outoforderRateList = new double[5];
public static double[] packetDelayList = new double[5];
public static double[] packetLossList = new double[5];


public static void main(String argv[]) {
	TotalDelay = 0;
	// RouterArrivalRate = Double.parseDouble(argv[1]);
	RouterServiceRate = 1250;
	SourceServiceRate = 1250;
	SourceServiceTime = 1/SourceServiceRate;
	NormalMean = Double.parseDouble(argv[0]);
	NormalStandarDeviation = Double.parseDouble(argv[1]);

	MeanInterArrivalTime = 1.0/1125; MeanServiceTime = 1.0/RouterServiceRate;
	SIGMA                = 0.6; TotalCustomers  = 300000;
	//long seed            = Long.parseLong(argv[0]);

	long seedList[] = {1,10,100,1000,10000};
	for (int k = 0; k<5; k++) {
	    //long seed            = Long.parseLong(argv[0]);
	    long seed = seedList[k];
		stream = new Random(seed);           // initialize rng stream
		
		RouterDelayTime = 0.05;
		SourceDelayTime = normal(stream, NormalMean, NormalStandarDeviation);
		FutureEventList = new EventList();
		HighQueue = new Queue();
		LowQueue = new Queue();
		packages = new ArrayList<Event>();
		Customers = new Queue();
		int dropped = 0;
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
		int currSequence = -1;
		Collections.sort(packages);
		int unOrderedCount = 0;

		for(int i=0; i<packages.size(); i++){
			// System.out.println("correct order is "+packages.get(i).get_order());
			/*
			if(i!=packages.get(i).get_order()){
				unOrderedCount++;
			}
			*/
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
		
		System.out.println("Before entering priority queue, out ot order rate is "+unOrderedCount/(double)TotalCustomers);
		
		System.out.println("");
		System.out.println("");
		System.out.println("");

		Clock = packages.get(0).get_arrive_router_time();
		LastEventTime = Clock;
		count = 0;
		NumberInService = 0;
		currSequence = -1;
		NumberOfDepartures = 0;
		FutureEventList = new EventList();
		Event e = new Event(arrival, packages.get(0).get_arrive_router_time(), count);
		// HighQueue.enqueue(e);
		FutureEventList.enqueue(e);
		// enqueueCount++;
		LowQueueLength = 0;
		HighQueueLength = 0;
		
		int j = 0;
		ArrayList<Event> arr = new ArrayList<Event>();
		finalArrivals = new ArrayList<Event>();
	  
	  	// Start to process the arrivals and departures on router
	  	while(NumberOfDepartures < TotalCustomers ) {
			// System.out.println("count "+count);
			// System.out.println("departures "+NumberOfDepartures);
			// System.out.println("enqueue count is "+enqueueCount);
			// System.out.println("getmin count is "+getMinCount);
			// System.out.println("");
		    Event evt = (Event)FutureEventList.getMin();  // get imminent event
			// getMinCount++;
		    FutureEventList.dequeue();                    // be rid of it
		    Clock = evt.get_time();                       // advance simulation time
			// System.out.println("event is "+evt.get_type()+" "+evt.get_order());
			// System.out.println("time is "+evt.get_time());
			// System.out.println("");
		
			if(j==0){
				currSequence = evt.get_order();
				arr.add(evt);
				ProcessRouterArrival(evt, 1);
				j++;
				continue;
			}

		    if( evt.get_type() == arrival ) {
				arr.add(evt);
				// System.out.println("order "+evt.get_order());
				if(evt.get_order()>currSequence){
					// System.out.println("cuurent sequence is "+currSequence);
					// System.out.println("low");
					if(LowQueueLength>=10000){
						
						dropped++;
						j++;
					    count++;
						if(count<TotalCustomers){
							Event next_arrival = new Event(arrival, packages.get(count).get_arrive_router_time(), packages.get(count).get_order());
							FutureEventList.enqueue( next_arrival );
							// enqueueCount++;
							LastEventTime = Clock;
						}
						
						continue;
					}else{
						currSequence = evt.get_order();
						if(HighQueueLength==0&&LowQueueLength==0&&NumberInService==0){
							ProcessRouterArrival(evt, 0);
						}else{
							LowQueue.enqueue(evt);
			//				System.out.println("lowqueue added "+evt.get_order());
			//				System.out.println("High order length "+HighQueueLength);
			//				System.out.println("Low order length "+LowQueueLength);
							LowQueueLength++;
							ProcessRouterArrival(evt, 0);
						}
					}
					
				} else{
			//		System.out.println("cuurent sequence is "+currSequence);
			//		System.out.println("high");
					if(HighQueueLength>=10000){
						
						dropped++;
						j++;
					    count++;
						if(count<TotalCustomers){
							Event next_arrival = new Event(arrival, packages.get(count).get_arrive_router_time(), packages.get(count).get_order());
							FutureEventList.enqueue( next_arrival );
					//		enqueueCount++;
							LastEventTime = Clock;
						}
						
						continue;
					} else{
						if(HighQueueLength==0&&LowQueueLength==0&&NumberInService==0){
							ProcessRouterArrival(evt, 1);
						}else{
							HighQueue.enqueue(evt);
				//			System.out.println("HHHHHHHHHHHHHHHHHHighqueue added "+evt.get_order());
			//				System.out.println("High order length "+HighQueueLength);
			//				System.out.println("Low order length "+LowQueueLength);
							HighQueueLength++;
							ProcessRouterArrival(evt, 1);
						}
						
					}
				}
			}

		    else  ProcessRouterDeparture(evt);
			j++;
	    }
		
		unOrderedCount = 0;
		//Collections.sort(finalArrivals);

		for(int i=0; i<finalArrivals.size(); i++){
			// System.out.println("chenged order "+finalArrivals.get(i).get_order());
			/*if(i!=finalArrivals.get(i).get_order()){
				unOrderedCount++;
			}
			if(arr.get(i).get_order()!=packages.get(i).get_order()){
				System.out.println("packages "+packages.get(i).get_order());
				System.out.println("arrival "+packages.get(i).get_order());
			}*/
			
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


		/********* newly added **********/
    	double outoforderRate = unOrderedCount/(double)TotalCustomers;
    	double packetDelay = TotalDelay/(double)TotalCustomers+0.05;
    	double packetLoss = dropped/(double)TotalCustomers;

    	outoforderRateList[k] = outoforderRate;
    	packetDelayList[k] = packetDelay;
    	packetLossList[k] = packetLoss;

	  	ReportGeneration();
	}


	double meanOutoforderRate = mean(outoforderRateList);
    double HOutoforderRate = ese(outoforderRateList);
    double CILowerOutoforderRate = meanOutoforderRate - HOutoforderRate;
    double CIUpperOutoforderRate = meanOutoforderRate + HOutoforderRate;
    
    System.out.println("CI lower end: \t" + CILowerOutoforderRate);
    System.out.println("mean: \t\t" + meanOutoforderRate);
    System.out.println("CI upper end: \t" + CIUpperOutoforderRate);
}



public static double mean(double numArray[]) {
	double sum = 0.0;
	int length = numArray.length;

	for(double num : numArray) {
	  sum += num;
	}

	double mean = sum/length;
	return mean;
}

// compute H = ese(sample mean)
public static double ese(double numArray[]) {
	double standardDeviation = 0.0;
	int length = numArray.length;
	double mean = mean(numArray);

	for(double num: numArray) {
	  standardDeviation += Math.pow(num - mean, 2);
	}
	standardDeviation = Math.sqrt(standardDeviation/(length-1));

	// t = 2.57, R = 5
	return (2.57 * standardDeviation / Math.sqrt(5));

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
  public static void ProcessRouterArrival(Event evt, int level) {

//  Customers.enqueue(evt); 
 // QueueLength++;
  // if the server is idle, fetch the event, do statistics
  // and put into service
  if( NumberInService == 0) {
	  if(LowQueueLength==0&&HighQueueLength==0){
	//	  System.out.println("aaaaaaaaaaasaa");
		//  finished = evt;
		  ScheduleRouterDeparture(evt);
	  }else{
		  ScheduleRouterDeparture(level);
	  }
	  
  }
  else TotalBusy += (Clock - LastEventTime);  // server is busy

  // adjust max queue length statistics
  if(MaxQueueLength < Math.max(HighQueueLength, LowQueueLength)) MaxQueueLength = Math.max(HighQueueLength, LowQueueLength);

  // schedule the next arrival
//  System.out.println("count is "+count);
  count++;
  if(count<TotalCustomers){
	  Event next_arrival = new Event(arrival, packages.get(count).get_arrive_router_time(), packages.get(count).get_order());
	  FutureEventList.enqueue( next_arrival );
//	  enqueueCount++;
	  LastEventTime = Clock; 
  }
 
 }

  public static void ScheduleRouterDeparture(Event evt) {
  double ServiceTime = MeanServiceTime;
  // get the job at the head of the queue
 // while (( ServiceTime = normal(stream, MeanServiceTime, SIGMA)) < 0 );
  //Event from low queue
 // ServiceTime = exponential(stream, MeanServiceTime);
  while (( ServiceTime = 1.0/1250.0) < 0 );
  Event depart = new Event(departure,Clock+ServiceTime, evt.get_order());
  finished = depart;
  totalWait+=(finished.get_time()-evt.get_time()-1.0/1250.0);
  TotalDelay += (finished.get_time()-evt.get_time());
  FutureEventList.enqueue( depart );
 // enqueueCount++;
  NumberInService = 1;
 }
 
 public static void ScheduleRouterDeparture(int level) {
  double ServiceTime = MeanServiceTime;
  // get the job at the head of the queue
 // while (( ServiceTime = normal(stream, MeanServiceTime, SIGMA)) < 0 );
  //Event from low queue
  while (( ServiceTime = 1.0/1250.0) < 0 );
  Event temp;
  if(HighQueueLength>0){
	  temp = (Event)HighQueue.dequeue();
	  //System.out.println("high!!!!!!!!!!!!!!!!!!!");
	  HighQueueLength--;
  } else {
	  temp = (Event)LowQueue.dequeue();
	  LowQueueLength--;
	//  System.out.println("LOWWWWWWWWWWWWWWWWWWWW");
  }
 // System.out.println("it's "+temp.get_order());
 // System.out.println("it's "+temp.get_time());
  Event depart = new Event(departure,Clock+ServiceTime, temp.get_order());
  finished = depart;
  totalWait+=(finished.get_time()-temp.get_time()-1.0/1250.0);
  TotalDelay += (finished.get_time()-temp.get_time());
 // System.out.println("its time "+depart.get_time());
  FutureEventList.enqueue( depart );
 // enqueueCount++;
  NumberInService = 1;
 }

public static void ProcessRouterDeparture(Event e) {
// System.out.println("departure "+e.get_order());
 // System.out.println("de time "+e.get_time());  
 // System.out.println("finished "+finished.get_order());	
 // System.out.println("fi time "+finished.get_time());
 // get the customer description
 //Event finished = (Event) Customers.dequeue();
 double delay = 0.05;
 TotalDelay += delay;
 //System.out.println("clock time is "+Clock);
 //System.out.println("delay is "+delay);
 //System.out.println("total time is "+(Clock+delay));
 finished.set_time(Clock+delay);
 finalArrivals.add(finished);
 // if there are customers in the queue then schedule
 // the departure of the next one
  if( HighQueueLength > 0 || LowQueueLength > 0 ) ScheduleRouterDeparture(1);
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
  //if (MaxQueueLength < QueueLength) MaxQueueLength = QueueLength;

  // schedule the next arrival
  Event next_arrival = new Event(arrival, Clock+exponential(stream, MeanInterArrivalTime), count);
  count++;
  FutureEventList.enqueue( next_arrival );
  LastEventTime = Clock;
 }

 public static void ScheduleDeparture() {
  double ServiceTime = MeanServiceTime;
  // get the job at the head of the queue
 // while (( ServiceTime = normal(stream, MeanServiceTime, SIGMA)) < 0 );
  while (( ServiceTime = 1.0/1250.0) < 0 );
  Event depart = new Event(departure,Clock+ServiceTime, count);
  FutureEventList.enqueue( depart );
  NumberInService = 1;
  QueueLength--;
 }

public static void ProcessDeparture(Event e) {
 // get the customer description
 Event finished = (Event) Customers.dequeue();
 double delay = normal(stream, NormalMean, NormalStandarDeviation);
 
 //System.out.println("clock time is "+Clock);
 //System.out.println("delay is "+delay);
 //System.out.println("total time is "+(Clock+delay));
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

