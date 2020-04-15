use std::str;
//use std::ascii::AsciiExt;
use std::time::{Duration, Instant};
use std::collections::HashMap;

// This file contains two simple structs for performance 
// timing measurements: TImedTask and PerfProfiler.
// PerfProfiler is a bit more simple to use, but has more
// overhead when starting and stopping timing events

// NOTE:
// Replace "TimedTask" with "DisabledTimedTask" to disable measurements


//TimedTask USAGE:
/*
    // Create and name regions in the code
    let mut t_chr =     TimedTask::new( "Loading ref chromosomes", "");
    let mut t_bam =     TimedTask::new( "Reading bam files", "");
    let mut t_pileups = TimedTask::new( "Processing pileups", "");
    //let mut t_build_indels = TimedTask::new( "Building indels", "Processing pileups");
    //let mut t_print_indels = TimedTask::new( "Printing indels", "Processing pileups");
    //let mut t_print_bs =     TimedTask::new( "Print base substitutions", "Processing pileups");
    let mut t_drain =   TimedTask::new( "Draining window", "");
    let mut t_total =   TimedTask::new( "Test total", "");
    let mut t_total2 =  TimedTask::new( "Test total2", "");
    let mut t_total3 =  TimedTask::new( "Test total3", "");

    let start_time = Instant::now(); //save processing start time
    t_total.start();
    t_total2.start();
    t_total3.start();

    // HEAVY PROCESSING

    t_total.end();
    t_total2.end();
    t_total3.end();
    eprintln!("INFO: Processing took {} seconds.", start_time.elapsed().as_secs());

    let mut task_list : Vec<TimedTask> = Vec::new();
    task_list.push( t_chr);
    task_list.push( t_bam);      
    task_list.push( t_pileups);
    task_list.push( t_drain);
    task_list.push( t_total);
    task_list.push( t_total2);
    task_list.push( t_total3);

    TimedTask::report( start_time, &mut task_list);
*/



pub struct TimedTask {
    name : String,
    dur : Duration,
    subtask_of : String,
    start : Instant,
    running : bool,
    startups : u64,
    stops : u64,    
}


pub struct DisabledTimedTask {
    running : bool,
}

impl DisabledTimedTask {

    pub fn new( _task_name : &str, _task_subtask_of : &str ) -> DisabledTimedTask {     
        DisabledTimedTask { running : false }        
    }
    pub fn is_enabled() -> bool { return false; }

    pub fn start(&mut self) -> () {}
    pub fn end(&mut self) -> () {}    
    fn drop(&mut self) {}
    fn add_duration( &mut self) -> () {}

    pub fn report( start_time : Instant, task_list : &mut Vec<DisabledTimedTask> ) -> () {
        return;
    }

}

impl TimedTask {

    pub fn new( task_name : &str, task_subtask_of : &str ) -> TimedTask {
     
        TimedTask {
            name : String::from( task_name),
            dur : Duration::new( 0u64, 0u32),
            subtask_of : String::from( task_subtask_of),
            start : Instant::now(),
            running : false,
            startups : 0u64,
            stops : 0u64,
        }        
    }

    pub fn is_enabled() -> bool { return true; }

    pub fn start(&mut self) -> () {

        //if self.running { self.add_duration(); }        
        if self.running { return; }
        self.startups += 1;
        self.running = true;
        self.start = Instant::now();
        
    }

    pub fn end(&mut self) -> () {

        if !self.running { return; }
        self.add_duration();
        self.running = false;
        self.stops += 1;
    }    

    fn drop(&mut self){
        //eprintln!( "Dropping TimedTask");
        self.end();
    }

    fn add_duration( &mut self) -> () {

        if !self.running { return; }
        self.dur += self.start.elapsed();               
    }


    pub fn report( start_time : Instant, task_list : &mut Vec<TimedTask> ) -> () {

        for tt in task_list.iter_mut() { tt.end(); }

        //use ansi_term::Colour::Green;  
        eprintln!( "\n-- TaskProfiler Report");
        eprintln!( "Measured task times:    (starts/stops)\n"); 

        let total_time = start_time.elapsed();
        eprintln!( "Total duration: {} seconds", duration_seconds( &total_time));

        //Print main tasks with subtasks
        for task in task_list.iter() {

                if !task.subtask_of.is_empty() { continue; } //Print non-subtasks
                TimedTask::print_task( &task, &task_list, total_time, 0usize);
                //eprintln!("{:5.1}% : {} ({:.2} seconds)", duration_percentage( &total_time, &task.dur), task.name, duration_seconds( &task.dur));            
        }
    }

    fn print_task( task : &TimedTask, task_list : &Vec<TimedTask>, total_time : Duration, level : usize) -> () {

        let mut level_str = String::new();

        // Create indentations for sub-tasks
        for _i in 0..level {
            level_str.push_str( "       > ");
            //level_str.push_str( "> ");
        }

        eprintln!("{}{:5.1}% : {} ({:.2} seconds) {}/{}", &level_str, duration_percentage( &total_time, &task.dur), task.name, duration_seconds( &task.dur), task.startups, task.stops);

        //Print subtasks recursively
        for t in task_list.iter() {

            //eprintln!( "sub:{}", t.subtask_of);
            if t.subtask_of == task.name {
                TimedTask::print_task( &t, task_list, task.dur, level+1 );
            }
            
        }
    }

}


pub struct PerfProfiler {
    
    start: Instant,
    tasks: HashMap<String, Duration>,
    cur_task : String,
    cur_time : Instant,

    sub_map : HashMap<String,Vec<String>>,
    sub_task : String,
    sub_time : Instant,

    enabled : bool,
}

fn duration_percentage( total : &Duration, dur : &Duration ) -> f32 {

    let s_total = (total.as_secs() as f64) + ((total.subsec_nanos() as f64) / 1000000000.0f64);
    let s_dur = (dur.as_secs() as f64) + ((dur.subsec_nanos() as f64) / 1000000000.0f64);
    return (s_dur / s_total * 100.0) as f32;
}

fn duration_seconds( dur : &Duration ) -> f32 {

    let s_dur = (dur.as_secs() as f64) + ((dur.subsec_nanos() as f64) / 1000000000.0f64);
    return s_dur as f32;
}


impl PerfProfiler {

    pub fn new() -> PerfProfiler {

        eprintln!("INFO: Measuring performance...");        

        PerfProfiler {
            start : Instant::now(),
            tasks : HashMap::new(),        
            cur_task : String::new(),
            cur_time : Instant::now(),
            sub_map : HashMap::new(),
            sub_task : String::new(),
            sub_time : Instant::now(),     
            enabled : true,       
        }        
    }

    fn drop(&mut self){
        self.report();
    }

    pub fn disable( &mut self) {
        self.enabled = false;
    }

    pub fn end_subtask( &mut self) -> () {

        if !self.enabled || self.sub_task.is_empty() { return; }        
        if let Some( x) = self.tasks.get_mut( &self.sub_task) {
            *x += self.sub_time.elapsed();        
        }
    }

    pub fn end_task( &mut self) -> () {

        if !self.enabled || self.cur_task.is_empty() { return; }
       
        //let task_opt = self.tasks.get_mut( &self.cur_task); //.or_insert( Duration::new());

        if let Some( x) = self.tasks.get_mut( &self.cur_task) {
            *x += self.cur_time.elapsed();
        }
        self.cur_task = String::new();
    }

    pub fn set_subtask( &mut self, task_name : &str ) -> () {

        if !self.enabled || task_name == &self.sub_task { return; }
  
        if self.cur_task.is_empty() { panic!("ERROR: main task not set for subtask {}.", &task_name); }

        self.end_subtask();     

        if task_name.is_empty() { return; }

        if !self.tasks.contains_key( task_name) {            
            
            //Add sub_task to full tasks list
            self.tasks.insert( String::from( task_name), Duration::new( 0u64, 0u32));            

            let mut pushed = false;
            if let Some( x) = self.sub_map.get_mut( &self.cur_task) {            
                x.push( String::from( task_name));
                pushed = true;
                //eprintln!( "Pushing: {}", &task_name);
            }

            if !pushed {
                let s = self.cur_task.clone(); //String::from( self.cur_task);
                let sub_vec = vec![String::from( task_name)];
                self.sub_map.insert( s, sub_vec);
                //eprintln!( "Inserting: {} -> {}", &self.cur_task, &task_name);
            }

        }

        self.sub_time = Instant::now();
        self.sub_task = String::from( task_name);


    }


	pub fn set_task( &mut self, task_name : &str ) -> () {

        if !self.enabled || task_name == &self.cur_task { return; }

        self.end_task();   
        self.end_subtask();     

        if task_name.is_empty() { return; }

        if !self.tasks.contains_key( task_name) {
            //self.tasks.insert( task_name, Duration:new( 0, 0));
            self.tasks.insert( String::from( task_name), Duration::new( 0u64, 0u32));
        }

        self.cur_time = Instant::now();
        self.cur_task = String::from( task_name);
    }    

    pub fn report( &mut self) -> () {

        if !self.enabled { return; }
        self.end_task();
        //use ansi_term::Colour::Red;
        //use ansi_term::Colour::Green;  
        eprintln!( "\nPerformance Profiler Report");
        eprintln!( "Measured task times:\n"); 

        let total_time = self.start.elapsed();
        eprintln!( "Total duration: {} seconds", duration_seconds( &total_time));

        let mut printed: HashMap<String,bool> = HashMap::new();

        //Print main tasks with subtasks
        for (task, vec) in &self.sub_map {

            if let Some( mt) = self.tasks.get( task) {
                let main_task_duration = duration_seconds( &mt);    
                eprintln!("{:5.1}% : {} ({:.2} seconds)", duration_percentage( &total_time, &mt), task, &main_task_duration);                                
                printed.insert( task.clone(), true);

                for sub_task_name in vec.iter() {
                    if let Some( st) = self.tasks.get( sub_task_name) {
                        let sub_task_duration = duration_seconds( &st);    
                        //let sub_dur = self.tasks.get( sub_task_name);
                        eprintln!("       > {:5.1}% : {} ({:.2} seconds)", duration_percentage( &mt, &st), sub_task_name, sub_task_duration);
                        printed.insert( sub_task_name.clone(), true);
                    }                
                }
            }
        }

        //Print main tasks that do not have subtasks
        for (task, duration) in &self.tasks {
            if printed.contains_key( task) { continue; }
            eprintln!("{:5.1}% : {} ({:.2} seconds)", duration_percentage( &total_time, &duration), task, duration_seconds( &duration));
        }

    }

}


