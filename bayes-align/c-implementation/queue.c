// a minimal implementation of queues via circular buffers for multithreading

// TODO: add mutexes to the postion indexes.

// FIFO queue via a ring buffer
typedef struct sb_queue
{
    int in;
    int out;
    int capacity;
    int length;
    sb_read queue[];
} sb_queue;


sb_queue new_sb_queue(int size)
{
    sb_queue.queue = sb_read array[size]
    sb_queue new_queue;

    new_queue.in = 0;
    new_queue.out = 0;
    new_queue.capacity = size;
    new_queue.length = 0;
    new_queue.queue = beginning;

    return(new_queue);
}

// add a read to an sb_queue, and returns 1 if there is an error
int enqueue(sb_queue *q, sb_read r)
{
    int position = q.in % q.capacity;

    // makes sure we are not about to clober the oldest living element in the buffer
    while (position == (q.out-1)%q.capacity)
    {
        continue;
    }

    q->queue[q.in % q.capacity] = sb_read;
    q.in = (q.in + 1) % q.capacity;
    q.length++;

}

// removes the longest living element in the queue
sb_read dequeue(sb_queue* q)
{

    int position = q.out % q.capacity;

    sb_read new_read;

    // TODO this needs to handle channel closures
    if new_read == NULL
    {
        return NULL;
    }

    new_read = q->array[position];

    q.out = (q.out + 1) % q.capacity;
    q.length--;

    return(new_read);
}
