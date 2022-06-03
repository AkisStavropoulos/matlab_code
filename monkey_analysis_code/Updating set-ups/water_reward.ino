const int led1 = 2; 
int x = 0; 

void setup() {
  Serial.begin(115200); 
  pinMode(led1, OUTPUT); 
}

void loop() {
  while(Serial.available()>0) {
    int set = Serial.parseInt();
    int del = Serial.parseInt(); 

    if (Serial.read() == '\n') {
      if (set == 1) {

        digitalWrite(led1, HIGH); 
        delay(del);
        digitalWrite(led1, LOW); 
      }
      }   
    }
}
